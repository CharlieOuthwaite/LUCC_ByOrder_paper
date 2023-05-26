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
moddir2 <- "5_RunLUIClimateModels/"
outDir <- "7_Predictions/Model_Comparisons/"
oridir1 <- "C:/Users/charl/Dropbox (UCL)/POSTDOC - BIOTA/0. PROJECTS/6. INSECTS LU CC/LanduseClimateInsects/2_RunSimpleLUIModel/"
oridir2 <- "C:/Users/charl/Dropbox (UCL)/POSTDOC - BIOTA/0. PROJECTS/6. INSECTS LU CC/LanduseClimateInsects/6_RunLUClimateModels/"

if(!dir.exists(outDir)) dir.create(outDir)


# load libraries
library(StatisticalModels)
library(predictsFunctions)
library(ggplot2)
source('0_Functions.R')


# read in the predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesClimate_Data.rds"))


##%######################################################%##
#                                                          #
####            Hyp 1: land use effect only             ####
#                                                          #
##%######################################################%##

#### Models in Outhwaite et al 2022) ####


##### Models in this study ####

load(file = paste0(moddir1, "Abundance_landuse_model.rdata")) # am3.3
load(file = paste0(moddir1, "Richness_landuse_model.rdata")) # sm3.3

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

# convert to percentage difference from primary vegetation
Coleoptera <- as.matrix(Coleoptera[, 1:1000])
col.preds <- sweep(x = Coleoptera, MARGIN = 2, STATS = Coleoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = col.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

colres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
colres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
colres$Order <- "Coleoptera"

# convert to percentage difference from primary vegetation
Diptera <- as.matrix(Diptera[, 1:1000])
dip.preds <- sweep(x = Diptera, MARGIN = 2, STATS = Diptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = dip.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

dipres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
dipres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
dipres$Order <- "Diptera"

# convert to percentage difference from primary vegetation
Hemiptera <- as.matrix(Hemiptera[, 1:1000])
hem.preds <- sweep(x = Hemiptera, MARGIN = 2, STATS = Hemiptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hem.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hemres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hemres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hemres$Order <- "Hemiptera"

# convert to percentage difference from primary vegetation
Hymenoptera <- as.matrix(Hymenoptera[, 1:1000])
hym.preds <- sweep(x = Hymenoptera, MARGIN = 2, STATS = Hymenoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hym.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hymres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hymres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hymres$Order <- "Hymenoptera"

# convert to percentage difference from primary vegetation
Lepidoptera <- as.matrix(Lepidoptera[, 1:1000])
lep.preds <- sweep(x = Lepidoptera, MARGIN = 2, STATS = Lepidoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = lep.preds,MARGIN = 1, FUN = median))*100)-100
grp.upper <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.975))*100)-100
grp.lower <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.025))*100)-100

lepres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
lepres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
lepres$Order <- "Lepidoptera"

# put it back together
result.sr <- rbind(colres, dipres, hemres, hymres, lepres)
result.sr$metric <- "species richness"

# now for the abundance model  
result.ab <- PredictGLMERRandIter(model = am3.3$model, data = data_tab)

# backtransform
result.ab <- exp(result.ab)-0.01

# convert to dataframe
result.ab <- as.data.frame(result.ab)

# add in the LU info
result.ab$LUI <- data_tab$LUI

# add in order info
result.ab$Order <- data_tab$Order

# express as a percentage of primary
# break into Orders
Order<- paste0("",result.ab$Order)
list.result.ab <- split(result.ab,Order)
list2env(list.result.ab,globalenv())

# convert to percentage difference from primary vegetation
Coleoptera <- as.matrix(Coleoptera[, 1:1000])
col.preds <- sweep(x = Coleoptera, MARGIN = 2, STATS = Coleoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = col.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

colres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
colres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
colres$Order <- "Coleoptera"

# convert to percentage difference from primary vegetation
Diptera <- as.matrix(Diptera[, 1:1000])
dip.preds <- sweep(x = Diptera, MARGIN = 2, STATS = Diptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = dip.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

dipres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
dipres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
dipres$Order <- "Diptera"

# convert to percentage difference from primary vegetation
Hemiptera <- as.matrix(Hemiptera[, 1:1000])
hem.preds <- sweep(x = Hemiptera, MARGIN = 2, STATS = Hemiptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hem.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hemres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hemres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hemres$Order <- "Hemiptera"

# convert to percentage difference from primary vegetation
Hymenoptera <- as.matrix(Hymenoptera[, 1:1000])
hym.preds <- sweep(x = Hymenoptera, MARGIN = 2, STATS = Hymenoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hym.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hymres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hymres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hymres$Order <- "Hymenoptera"

# convert to percentage difference from primary vegetation
Lepidoptera <- as.matrix(Lepidoptera[, 1:1000])
lep.preds <- sweep(x = Lepidoptera, MARGIN = 2, STATS = Lepidoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = lep.preds,MARGIN = 1, FUN = median))*100)-100
grp.upper <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.975))*100)-100
grp.lower <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.025))*100)-100

lepres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
lepres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
lepres$Order <- "Lepidoptera"

# put it back together
result.ab <- rbind(colres, dipres, hemres, hymres, lepres)
result.ab$metric <- "total abundance"

all_res <- rbind(result.ab, result.sr)

# save table
write.csv(all_res, file = paste0(outDir, "/percentage_change_LU_Order.csv"))


##%######################################################%##
#                                                          #
####        Figure for percentage changes by LU         ####
#                                                          #
##%######################################################%##

# load version including results from Outhwaite et al 2022 and rounded values to 2dp
all_res <- read.csv(file = paste0(outDir, "/percentage_change_LU_Order_inc_2022.csv"))

all_res$LUI <- factor(all_res$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Low-intensity agriculture", "High-intensity agriculture"))
all_res$Order <- factor(all_res$Order, levels = c("All insects", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"))


# create point and error bar plot
ggplot(data = all_res, aes(col = LUI, group = Order)) + 
  geom_point(aes(x = LUI, y = perc, shape = Order), size = 1.5, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = LUI, ymin = lower_CI, ymax = upper_CI), position= position_dodge(width = 1), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  facet_wrap(~ Metric) +
  xlab("") +
  scale_y_continuous(limits = c(-80, 80), breaks = scales::pretty_breaks(n = 10)) +
  ylab("Percentage change (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00"), guide = "none") +
  scale_shape_manual(values=c(16, 17, 18, 15, 0, 1))+
    theme(legend.position = "bottom", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2), 
        text = element_text(size = 7), 
        legend.key=element_blank(), 
        legend.title = element_blank())


ggsave(filename = paste0(outDir, "Comparison_LU_only.pdf"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "Comparison_LU_only.jpeg"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)




##%######################################################%##
#                                                          #
####  Hyp 2: Land use and climate anomaly interaction   ####
#                                                          #
##%######################################################%##

#### Outhwaite et al 2022 ####

predictsSites <- readRDS(paste0(oridir2,"PREDICTSSiteData.rds"))

# load in models
load(file = paste0(oridir2, "MeanAnomalyModelAbund.rdata")) # MeanAnomalyModelAbund
load(file = paste0(oridir2, "MeanAnomalyModelRich.rdata")) # MeanAnomlayModelRich


# what is the rescaled value of STA of 1
BackTransformCentreredPredictor(transformedX = 0.97, originalX = predictsSites$StdTmeanAnomaly) # 0.97 gives about 1 

# what is the rescaled value of STA of 0
BackTransformCentreredPredictor(transformedX = -1.39, originalX = predictsSites$StdTmeanAnomaly) # -1.39 gives about 0 

# set up table for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.97),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
# row 1

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to dataframe
#result.ab <- as.data.frame(a.preds.tmean)

# # add in the LU info
# result.ab$LUI <- nd$UI2

# convert to percentage difference from primary vegetation
a.preds <- sweep(x = a.preds.tmean, MARGIN = 2, STATS = a.preds.tmean[1,], FUN = '/')

# get quantiles
a.preds.median <- ((apply(X = a.preds,MARGIN = 1,FUN = median))*100)-100
a.preds.upper <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
a.preds.lower <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


## species richness predictions ##
s.preds <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model, data = nd)

s.preds <- exp(s.preds)

# convert to percentage difference from primary vegetation
s.preds <- sweep(x = s.preds, MARGIN = 2, STATS = s.preds[1,], FUN = '/')

# get quantiles
s.preds.median <- ((apply(X = s.preds,MARGIN = 1,FUN = median))*100)-100
s.preds.upper <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
s.preds.lower <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


# combine data into one table for plotting
abun_res <- as.data.frame(cbind(a.preds.median, a.preds.lower, a.preds.upper))
rich_res <- as.data.frame(cbind(s.preds.median, s.preds.lower, s.preds.upper))
colnames(abun_res) <- c("perc", "lower_CI", "upper_CI")
colnames(rich_res) <- c("perc", "lower_CI", "upper_CI")
abun_res$Metric <- "Total abundance"
rich_res$Metric <- "Species richness"
abun_res$LUI <- nd$UI2
rich_res$LUI <- nd$UI2

all_res_ori <- rbind(abun_res, rich_res)
all_res_ori$Study <- "Outhwaite et al 2022"

all_res_ori$Fixed_effs <- "Land use and climate"

all_res_ori$Order <- "All insects"

all_res_ori$Realm <- "Global"

all_res_ori$STA <- c(0,1)

all_res_ori <- all_res_ori[, c(6, 4, 8, 7, 9, 5, 1:3, 10)]


#### This study ####

# need predictions for each land use for 0 and 1 STA for each order

# read in the predicts data

# load in models
load(file = paste0(moddir2, "MeanAnomalyModelAbund.rdata")) # MeanAnomalyModelAbund
load(file = paste0(moddir2, "MeanAnomalyModelRich.rdata")) # MeanAnomlayModelRich

#### abundance predictions ####

# what is the rescaled value of STA of 1
BackTransformCentreredPredictor(transformedX = 0.999, originalX = predictsSites$StdTmeanAnomaly) # 0.999 gives about 1 

# what is the rescaled value of STA of 0
BackTransformCentreredPredictor(transformedX = -1.39, originalX = predictsSites$StdTmeanAnomaly) # -1.39 gives about 0 

# set up table for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.999),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
# row for each order

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to dataframe
result.ab <- as.data.frame(a.preds.tmean)

# add in the LU info
result.ab$LUI <- nd$LUI

# add in order info
result.ab$Order <- nd$Order

# express as a percentage of primary
# break into Orders
Order<- paste0("",result.ab$Order)
list.result.ab <- split(result.ab,Order)
list2env(list.result.ab,globalenv())

# convert to percentage difference from primary vegetation
Coleoptera <- as.matrix(Coleoptera[, 1:1000])
col.preds <- sweep(x = Coleoptera, MARGIN = 2, STATS = Coleoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = col.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

colres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
colres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
colres$Order <- "Coleoptera"

# convert to percentage difference from primary vegetation
Diptera <- as.matrix(Diptera[, 1:1000])
dip.preds <- sweep(x = Diptera, MARGIN = 2, STATS = Diptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = dip.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

dipres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
dipres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
dipres$Order <- "Diptera"

# convert to percentage difference from primary vegetation
Hemiptera <- as.matrix(Hemiptera[, 1:1000])
hem.preds <- sweep(x = Hemiptera, MARGIN = 2, STATS = Hemiptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hem.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hemres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hemres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
hemres$Order <- "Hemiptera"

# convert to percentage difference from primary vegetation
Hymenoptera <- as.matrix(Hymenoptera[, 1:1000])
hym.preds <- sweep(x = Hymenoptera, MARGIN = 2, STATS = Hymenoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hym.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hymres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hymres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
hymres$Order <- "Hymenoptera"

# convert to percentage difference from primary vegetation
Lepidoptera <- as.matrix(Lepidoptera[, 1:1000])
lep.preds <- sweep(x = Lepidoptera, MARGIN = 2, STATS = Lepidoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = lep.preds,MARGIN = 1, FUN = median))*100)-100
grp.upper <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.975))*100)-100
grp.lower <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.025))*100)-100

lepres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
lepres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
lepres$Order <- "Lepidoptera"

# put it back together
result.ab <- rbind(colres, dipres, hemres, hymres, lepres)
result.ab$Metric <- "Total abundance"

#### species richness predictions ####

# set up table for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.999),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
# row for each order, its the 1st row each time

# predict the results
s.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd)

# back transform the abundance values
s.preds.tmean <- exp(s.preds.tmean)

# convert to dataframe
result.sr <- as.data.frame(s.preds.tmean)

# add in the LU info
result.sr$LUI <- nd$LUI

# add in order info
result.sr$Order <- nd$Order

# express as a percentage of primary
# break into Orders
Order<- paste0("",result.sr$Order)
list.result.sr <- split(result.sr,Order)
list2env(list.result.sr,globalenv())

# convert to percentage difference from primary vegetation
Coleoptera <- as.matrix(Coleoptera[, 1:1000])
col.preds <- sweep(x = Coleoptera, MARGIN = 2, STATS = Coleoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = col.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

colres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
colres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
colres$Order <- "Coleoptera"

# convert to percentage difference from primary vegetation
Diptera <- as.matrix(Diptera[, 1:1000])
dip.preds <- sweep(x = Diptera, MARGIN = 2, STATS = Diptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = dip.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

dipres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
dipres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
dipres$Order <- "Diptera"

# convert to percentage difference from primary vegetation
Hemiptera <- as.matrix(Hemiptera[, 1:1000])
hem.preds <- sweep(x = Hemiptera, MARGIN = 2, STATS = Hemiptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hem.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hemres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hemres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
hemres$Order <- "Hemiptera"

# convert to percentage difference from primary vegetation
Hymenoptera <- as.matrix(Hymenoptera[, 1:1000])
hym.preds <- sweep(x = Hymenoptera, MARGIN = 2, STATS = Hymenoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hym.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hymres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hymres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
hymres$Order <- "Hymenoptera"

# convert to percentage difference from primary vegetation
Lepidoptera <- as.matrix(Lepidoptera[, 1:1000])
lep.preds <- sweep(x = Lepidoptera, MARGIN = 2, STATS = Lepidoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = lep.preds,MARGIN = 1, FUN = median))*100)-100
grp.upper <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.975))*100)-100
grp.lower <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.025))*100)-100

lepres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
lepres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
lepres$Order <- "Lepidoptera"

# put it back together
result.sr <- rbind(colres, dipres, hemres, hymres, lepres)
result.sr$Metric <- "Species Richness"


all_res <- rbind(result.ab, result.sr)

all_res$Study <- "This study"

all_res$Fixed_effs <- "Land use and climate"

all_res$Realm <- "Global"

all_res$STA <- c(0, 1)

all_res <- all_res[ , c(7,6,5, 8, 9, 4, 1:3, 10)]

names(all_res)[7:9] <- c("perc", "lower_CI", "upper_CI")

final_res <- rbind(all_res_ori, all_res)

# save table
write.csv(final_res, file = paste0(outDir, "/percentage_change_LU_STA_Order_inc2022.csv"))


##%######################################################%##
#                                                          #
####     Figure for percentage changes by LU and STA    ####
#                                                          #
##%######################################################%##

# load version including results from Outhwaite et al 2022 and rounded values to 2dp
final_res <- read.csv(file = paste0(outDir, "/percentage_change_LU_STA_Order_inc2022.csv"))

final_res$Order <- factor(final_res$Order, levels = c("All insects", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"))

final_res$STA <- factor(final_res$STA, levels = c("0", "1"))

final_res$LUI <- factor(final_res$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

plot_data <- final_res[final_res$Metric == "Total abundance", ]


# create point and error bar plot
ggplot(data = plot_data, aes(col = LUI, group = STA)) + 
  geom_point(aes(x = LUI, y = perc, shape = STA), size = 1.5, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = LUI, ymin = lower_CI, ymax = upper_CI), position= position_dodge(width = 1), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  facet_wrap(~ Order) +
  xlab("") +
  scale_y_continuous(limits = c(-100, 120), breaks = scales::pretty_breaks(n = 10)) +
  ylab("Percentage change in total abundance (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00"), guide = "none") +
  scale_shape_manual(values=c(16, 17, 18, 15, 0, 1))+
  theme(legend.position = "bottom", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2), 
        text = element_text(size = 7), 
        legend.key=element_blank())


ggsave(filename = paste0(outDir, "Comparison_LUSTA_Abun.pdf"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "Comparison_LUSTA_Abun.jpeg"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)

### now richness ###
plot_data <- final_res[final_res$Metric == "Species richness" | final_res$Metric == "Species Richness",]

# create point and error bar plot
ggplot(data = plot_data, aes(col = LUI, group = STA)) + 
  geom_point(aes(x = LUI, y = perc, shape = STA), size = 1.5, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = LUI, ymin = lower_CI, ymax = upper_CI), position= position_dodge(width = 1), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  facet_wrap(~ Order) +
  xlab("") +
  scale_y_continuous(limits = c(-100, 270), breaks = scales::pretty_breaks(n = 10)) +
  ylab("Percentage change in species richness (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00"), guide = "none") +
  scale_shape_manual(values=c(16, 17, 18, 15, 0, 1))+
  theme(legend.position = "bottom", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2), 
        text = element_text(size = 7), 
        legend.key=element_blank())


ggsave(filename = paste0(outDir, "Comparison_LUSTA_Rich.pdf"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "Comparison_LUSTA_Rich.jpeg"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)

