############################################################
#                                                          #
#           Predictions from both sets of models           #
#                                                          #
############################################################


# In this script, I predict the biodiversity for different combinations of
# land use and STA for both the models from Outhwaite et al (2022) and those
# tested here. 

rm(list = ls())

# directories
predictsDataDir <- "5_RunLUIClimateModels/"
moddir1 <- "2_RunSimpleLUIModel/"
moddir2 <- "5_LUIClimateModels/"
outDir <- "7_Predictions/Model_Comparisons/"
oridir1 <- "C:/Users/charl/Dropbox (UCL)/POSTDOC - BIOTA/0. PROJECTS/6. INSECTS LU CC/LanduseClimateInsects/6_RunLUClimateModels"
oridir2 <- "C:/Users/charl/Dropbox (UCL)/POSTDOC - BIOTA/0. PROJECTS/6. INSECTS LU CC/LanduseClimateInsects/2_RunSimpleLUIModel"

if(!dir.exists(outDir)) dir.create(outDir)


# load libraries
library(StatisticalModels)
library(predictsFunctions)
source('0_Functions.R')


##### 1. This Study #####

# read in the predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesClimate_Data.rds"))


###### Hyp 1: land use effect only ######

# these values can be retrieved from the info in script 2_RunSimpleLUIModel.R
# the median values were also used for plotting of Figure 1.

load(file = paste0(moddir1, "Abundance_landuse_model.rdata"))
load(file = paste0(moddir1, "Richness_landuse_model.rdata"))

data_tab <- expand.grid(LUI = factor(c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"), levels = levels(am3.3$data$LUI)), 
                       Order = factor(c("Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"), levels = levels(am3.3$data$Order)),
                       LogAbund = 0,
                       Species_richness = 0)


# predict the results
result.sr <- PredictGLMERRandIter(model = sm3.3$model, data = data_tab)

# backtransform
result.sr <- exp(result.sr)

# convert to dataframe
result.sr <- as.data.frame(result.sr)

# add in the LUI info
result.sr$LUI <- data_tab$LUI

# add in the Order info
result.sr$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr$Order)
list.result.sr <- split(result.sr,Order)
list2env(list.result.sr,globalenv())

results_ab <- NULL

# convert to percentage difference from primary vegetation
Coleoptera <- as.matrix(Coleoptera[, 1:1000])
col.preds <- sweep(x = Coleoptera, MARGIN = 2, STATS = Coleoptera[1,], FUN = '/')

# get quantiles
col.median <- ((apply(X = col.preds,MARGIN = 1,FUN = median))*100)-100
col.upper <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
col.lower <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

colres <- as.data.frame(cbind(col.median, col.upper, col.lower))
colres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
colres$Order <- "Coleoptera"

#### UP TO HERE ####

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100

# put it back together
result.ab <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera)



# express as a percentage of primary
result.sr$perc <- ((result.sr$y/result.sr$y[1]) * 100) - 100

# add in SCA vals
result.sr$SCA <- c(0, 1, 1, 0,0)  

# now for the abundance model  
result.ab <- PredictGLMER(model = MeanAnomalyModelAbund$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab <- exp(result.ab)-0.01

# add in the LU info
result.ab$UI2 <- data_tab$UI2

# express as a percentage of primary
result.ab$perc <- ((result.ab$y/result.ab$y[1]) * 100) - 100

# add in SCA vals
result.ab$SCA <- c(0, 1, 1, 0,0) 

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)

all_res$measure <- c(rep("ab", 5), rep("sr", 5))

# save table
write.csv(all_res, file = paste0(outDir, "/percentage_change_LU_CC.csv"))


###### Hyp 2: Land use and climate anomaly interaction ######

# first one, looking at SCA of 1

# load in models
load(file = paste0(moddir, "/MeanAnomalyModelAbund.rdata"))
load(file = paste0(moddir, "/MeanAnomalyModelRich.rdata"))
#load(file = paste0(moddir, "/MaxAnomalyModelAbund.rdata"))
#load(file = paste0(moddir, "/MaxAnomalyModelRich.rdata"))

# create matrix for predictions
# Primary, Low, High
# SCA = 1
# abun and richness = 0

# what is the rescaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = 0.98, originalX = predictsSites$StdTmeanAnomaly) # 0.98 gives about 1 

# what is the rescaled value of SCA of 0
BackTransformCentreredPredictor(transformedX = -1.39, originalX = predictsSites$StdTmeanAnomaly) # -1.39 gives about 0 

# reference is primary with 0 climate change so have 0 for that row

data_tab <- data.frame(UI2 = c("Primary vegetation", "Agriculture_Low", "Agriculture_High", "Agriculture_Low", "Agriculture_High"), 
                       StdTmeanAnomalyRS = c(-1.39,0.98,0.98,-1.39,-1.39),
                       LogAbund = 0,
                       Species_richness = 0)

# factor the LU info
data_tab$UI2 <- factor(data_tab$UI2, levels = levels(predictsSites$UI2))

# predict the results
result.sr <- PredictGLMER(model = MeanAnomalyModelRich$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr <- exp(result.sr)

# add in the LU info
result.sr$UI2 <- data_tab$UI2

# express as a percentage of primary
result.sr$perc <- ((result.sr$y/result.sr$y[1]) * 100) - 100

# add in SCA vals
result.sr$SCA <- c(0, 1, 1, 0,0)  

# now for the abundance model  
result.ab <- PredictGLMER(model = MeanAnomalyModelAbund$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab <- exp(result.ab)-0.01

# add in the LU info
result.ab$UI2 <- data_tab$UI2

# express as a percentage of primary
result.ab$perc <- ((result.ab$y/result.ab$y[1]) * 100) - 100

# add in SCA vals
result.ab$SCA <- c(0, 1, 1, 0,0) 

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)

all_res$measure <- c(rep("ab", 5), rep("sr", 5))

# save table
write.csv(all_res, file = paste0(outDir, "/percentage_change_LU_CC.csv"))


### Hyp 3:  Land use, climate and natural habitat interactions ###

# 25% NH cover = -1.05 in rescaled values
# 50% NH cover = 0.01493712 in rescaled values
# 75% NH cover = 1.045653 in rescaled values
# 100% NH cover = 2.073849 in rescaled values

# what is the reference here
# primary vegetation, SCA = 0, NH = 100

# what is the resscaled value of NH
#BackTransformCentreredPredictor(transformedX = 2.073849, originalX = predictsSites$NH_5000) # -0.99 gives about 0 


load(file = "7_RunLUClimateNHModels/MeanAnomalyModelAbun_NH.rdata")
load(file = "7_RunLUClimateNHModels/RichMeanAnomalyModel_NH.rdata")



data_tab <- data.frame(UI2 = c("Primary vegetation", rep("Agriculture_Low", 2), rep("Agriculture_High", 2)), 
                       StdTmeanAnomalyRS = c(-0.95,1.66,1.66, 1.66,1.66),
                       NH_5000.rs = c(2.073849, 1.045653, -1.05, 1.045653, -1.05 ), #100, 75, 25
                       LogAbund = 0,
                       Species_richness = 0)

# factor the LU info
data_tab$UI2 <- factor(data_tab$UI2, levels = levels(predictsSites$UI2))


# now for the abundance model  
result.ab <- PredictGLMER(model = AbundMeanAnomalyModel1$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab <- exp(result.ab)-0.01

# add in the LU info
result.ab$UI2 <- data_tab$UI2

result.ab$NH <- c(100, 75, 25, 75, 25)
result.ab$CC <- c(0, 1, 1, 1, 1)

# express as a percentage of primary
result.ab$perc <- ((result.ab$y/result.ab$y[1]) * 100) - 100

result.ab$metric <- "ab"


result.sr <- PredictGLMER(model = RichMeanAnomalyModel1$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr <- exp(result.sr)

# add in the LU info
result.sr$UI2 <- data_tab$UI2

result.sr$NH <- c(100, 75, 25, 75, 25)
result.sr$CC <- c(0, 1, 1, 1, 1)

# express as a percentage of primary
result.sr$perc <- ((result.sr$y/result.sr$y[1]) * 100) - 100

result.sr$metric <- "sr"

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)

# save table
write.csv(all_res, file = paste0(outDir, "/percentage_change_LU_CC_NH.csv"))



