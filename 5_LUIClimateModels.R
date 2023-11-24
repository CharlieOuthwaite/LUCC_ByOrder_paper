##%######################################################%##
#                                                          #
####    Run models for climate land use interactions    ####
#                         Global                           #
#                                                          #
##%######################################################%##

# in this script, models are run to assess the response of insect 
# biodiversity to the land use intensity metric and to the climate
# anomalies and their interaction. 

# clear working environment
rm(list = ls())

# set directories 
inDir <- "4_PREDICTSMatchClimateIndex/"
outDir <- "5_RunLUIClimateModels/"
plotDir <- "5_RunLUIClimateModels/Plots/"
tabDir <- "5_RunLUIClimateModels/Tables/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(plotDir)) dir.create(plotDir)
if(!dir.exists(tabDir)) dir.create(tabDir)

# sink(paste0(outDir,"log_LUI_ClimateModels.txt"))
# 
# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
packages_model <- c("devtools","StatisticalModels", "predictsFunctions","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "ggplot2", "cowplot", "sjPlot", "yarg", "lme4", "gt", "broom.mixed", "MASS")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")

#### 1. Organise data ####

# read in the predicts data
predictsSites <- readRDS(paste0(inDir,"PREDICTSSites_Climate.rds"))
predictsSites <- predictsSites@data

# set LUI as factor and set reference level
predictsSites$LUI <- factor(predictsSites$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# rescale the climate variables
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)

predictsSites <- droplevels(predictsSites)

# some of the climate values are NA since they do not meet the thresholds
predictsSites <- predictsSites[!is.na(predictsSites$avg_temp), ] # 8858 rows

# take a look at possible correlations between variables
cor(predictsSites$avg_temp, predictsSites$TmeanAnomaly)

# -0.2657298

cor(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly)

# 0.1973725

cor(predictsSites$TmeanAnomaly, predictsSites$StdTmeanAnomaly)

# 0.2519302

# save the dataset
saveRDS(object = predictsSites, file = paste0(outDir, "PREDICTSSitesClimate_Data.rds"))
#predictsSites <- readRDS(file = paste0(outDir, "PREDICTSSitesClimate_Data.rds"))


#### 2. Create models of climate land use interactions ####

# i. Abundance, mean anomaly including interaction

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ] # 8468 rows

# run model
MeanAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS", 'Latitude', 'Longitude'))

# get summary
summary(MeanAnomalyModelAbund$model)

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "MeanAnomalyModelAbund.rdata"))


# i. Abundance, mean anomaly excluding interaction with order

MeanAnomalyModelAbund2 <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "LUI * StdTmeanAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS", 'Latitude', 'Longitude'))

# get summary
summary(MeanAnomalyModelAbund2$model)


# save the model output
save(MeanAnomalyModelAbund2, file = paste0(outDir, "MeanAnomalyModelAbund_noOrder.rdata"))


# ii. Richness, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 8858 rows

MeanAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS", 'Latitude', 'Longitude'))

# get summary
summary(MeanAnomalyModelRich$model)

# save the model output
save(MeanAnomalyModelRich, file = paste0(outDir, "MeanAnomalyModelRich.rdata"))

# ii. Richness, mean anomaly excluding interaction with order

MeanAnomalyModelRich2 <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "LUI * StdTmeanAnomalyRS",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS", 'Latitude', 'Longitude'))

# save the model output
summary(MeanAnomalyModelRich2$model)

# save the model output
save(MeanAnomalyModelRich2, file = paste0(outDir, "MeanAnomalyModelRich_noOrder.rdata"))


# save model output tables for use in supplementary information
  # use function from sjPlot library to save neat versions of model output table
  # conditional: the conditional R2 value, i.e. the variance explained by fixed and random effects 
  # marginal: the marginal R2 value, i.e. the variance explained by the fixed effects

tab_model(MeanAnomalyModelAbund$model, transform = NULL, file = paste0(tabDir,"AbunMeanAnom_output_table.html"))
summary(MeanAnomalyModelAbund$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund$model) # check the R2 values 
# $conditional
# [1] 0.4164462
# 
# $marginal
# [1] 0.08326023

tab_model(MeanAnomalyModelRich$model, transform = NULL, file = paste0(tabDir,"RichMeanAnom_output_table.html"))
summary(MeanAnomalyModelRich$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich$model) # check the R2 values
# $conditional
# [1] 0.7529901
# 
# $marginal
# [1] 0.1093307

tab_model(MeanAnomalyModelAbund2$model, transform = NULL, file = paste0(tabDir,"AbunMeanAnom_output_table_noOrder.html"))
summary(MeanAnomalyModelAbund2$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund2$model) # check the R2 values 

# $conditional
# [1] 0.3417518
# 
# $marginal
# [1] 0.03110523

tab_model(MeanAnomalyModelRich2$model, transform = NULL, file = paste0(tabDir,"RichMeanAnom_output_table_noOrder.html"))
summary(MeanAnomalyModelRich2$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich2$model) # check the R2 values

# $conditional
# [1] 0.6715988
# 
# $marginal
# [1] 0.009126409


##%######################################################%##
#                                                          #
####           Richness and abundance plots             ####
#                          Global                          #
#                                                          #
##%######################################################%##

# load dataset and models if not already in environment
# predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSitesClimate_Data.rds"))
# load(paste0(outDir, "MeanAnomalyModelAbund.rdata"))
# load(paste0(outDir, "MeanAnomalyModelRich.rdata"))

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

#### 3. Plots of results for abundance, mean Anomaly ####

nd <- expand.grid(
      StdTmeanAnomalyRS = seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        length.out = 100),
      LUI = factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(MeanAnomalyModelAbund$data$LUI)),
      Order = factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
                      transformedX = nd$StdTmeanAnomalyRS,
                      originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# set quantiles by Order
# coleoptera
C_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)

# diptera
D_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)

# hemiptera
He_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)

# hymenoptera
Hy_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)

# lepidoptera
L_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)


# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model, data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# create list of matrices
number_of_chunks = 5

list_a.preds.tmean <- lapply(seq(1, NROW(a.preds.tmean), ceiling(NROW(a.preds.tmean)/number_of_chunks)),
       function(i) a.preds.tmean[i:min(i + ceiling(NROW(a.preds.tmean)/number_of_chunks) - 1, NROW(a.preds.tmean)),])

# name them
names(list_a.preds.tmean) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")
                             
# sweep out according to the refRow
list_a.preds.tmean <- lapply(list_a.preds.tmean,FUN=function(x){
 sweep (x=x, MARGIN = 2, STATS=x[refRow[1],],FUN="/") 
})

# print to global environment
list2env(list_a.preds.tmean,globalenv())

# split nd by order
Order <- paste0("nd_",nd$Order)

# create a list of data frames
by_Order <- split(nd,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QPV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QPV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QSV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QSV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QAL[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QAL[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QAH[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QAH[2])),] <- NA

Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < D_QPV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > D_QPV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < D_QSV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > D_QSV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS < D_QAL[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS > D_QAL[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS < D_QAH[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS > D_QAH[2])),] <- NA

Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QPV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QPV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QSV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QSV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QAL[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QAL[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QAH[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QAH[2])),] <- NA

Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QPV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QPV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QSV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QSV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QAL[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QAL[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QAH[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QAH[2])),] <- NA

Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QPV[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QPV[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Secondary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QSV[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Secondary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QSV[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_Low") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QAL[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_Low") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QAL[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_High") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QAH[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_High") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot
nd_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = median,na.rm=TRUE))*100)-100
nd_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = median,na.rm=TRUE))*100)-100
nd_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                  FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                  FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd_Coleoptera$LUI <- factor(nd_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Diptera$LUI <- factor(nd_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hemiptera$LUI <- factor(nd_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hymenoptera$LUI <- factor(nd_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Lepidoptera$LUI <- factor(nd_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Coleoptera$PredLower, ymax = nd_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") + 
  ggtitle("Coleoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 
  

p_diptera <- ggplot(data = nd_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Diptera$PredLower, ymax = nd_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Diptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 
  

p_hemiptera <- ggplot(data = nd_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hemiptera$PredLower, ymax = nd_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hemiptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))
  

p_hymenoptera <- ggplot(data = nd_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hymenoptera$PredLower, ymax = nd_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hymenoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_lepidoptera <- ggplot(data = nd_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Lepidoptera$PredLower, ymax = nd_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Lepidoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 
  

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 7), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomAbund <- cowplot::plot_grid(p_coleoptera, p_diptera, p_hemiptera, p_hymenoptera, p_lepidoptera, legend, ncol=3)

# save plot (pdf)
ggsave(filename = paste0(plotDir, "MeanAnomAbund.pdf"), plot = MeanAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("MeanAnomAbund.jpeg", device ="jpeg", path = plotDir, width=20, height=15, units="cm", dpi = 350)


#### 4. Richness, Mean Anomaly ####

nd2 <- expand.grid(StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                                         to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                                         length.out = 100),
                   LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                              levels = levels(MeanAnomalyModelRich$data$LUI)),
                   Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$LUI=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# set quantiles
# coleoptera
C_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)

# diptera
D_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)

# hemiptera
He_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)

# hymenoptera
Hy_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)

# lepidoptera
L_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)


# predict the results
sr.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model, data = nd2)

# back transform the abundance values
sr.preds.tmean <- exp(sr.preds.tmean)

# create list of matrices
number_of_chunks = 5

list_sr.preds.tmean <- lapply(seq(1, NROW(sr.preds.tmean), ceiling(NROW(sr.preds.tmean)/number_of_chunks)),
                              function(i) sr.preds.tmean[i:min(i + ceiling(NROW(sr.preds.tmean)/number_of_chunks) - 1, NROW(sr.preds.tmean)),])
# name them
names(list_sr.preds.tmean) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")

# sweep out according to the refRow
list_sr.preds.tmean <- lapply(list_sr.preds.tmean,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[refRow[1],],FUN="/") 
})
                              
list2env(list_sr.preds.tmean,globalenv())

# split nd by order
Order<- paste0("nd2_",nd2$Order)

# create a list of data frames
by_Order <- split(nd2,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QPV[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QPV[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QSV[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QSV[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QAL[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QAL[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QAH[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QAH[2])),] <- NA

Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < D_QPV[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > D_QPV[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < D_QSV[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > D_QSV[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS < D_QAL[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS > D_QAL[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS < D_QAH[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS > D_QAH[2])),] <- NA

Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QPV[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QPV[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QSV[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QSV[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QAL[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QAL[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QAH[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QAH[2])),] <- NA

Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QPV[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QPV[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QSV[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QSV[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QAL[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QAL[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QAH[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QAH[2])),] <- NA

Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QPV[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QPV[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QSV[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QSV[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QAL[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QAL[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QAH[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QAH[2])),] <- NA



# Get the median, upper and lower quants for the plot
nd2_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd2_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                  FUN = median,na.rm=TRUE))*100)-100
nd2_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd2_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd2_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd2_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd2_Coleoptera$LUI <- factor(nd2_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Diptera$LUI <- factor(nd2_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Hemiptera$LUI <- factor(nd2_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Hymenoptera$LUI <- factor(nd2_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Lepidoptera$LUI <- factor(nd2_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot
p_coleoptera <- ggplot(data = nd2_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Coleoptera$PredLower, ymax = nd2_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Coleoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 

p_diptera <- ggplot(data = nd2_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Diptera$PredLower, ymax = nd2_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Diptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_hemiptera <- ggplot(data = nd2_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hemiptera$PredLower, ymax = nd2_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hemiptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_hymenoptera <- ggplot(data = nd2_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hymenoptera$PredLower, ymax = nd2_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hymenoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_lepidoptera <- ggplot(data = nd2_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Lepidoptera$PredLower, ymax = nd2_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Lepidoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 7), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomRich <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,legend, ncol=3)

# save plot (pdf)
ggsave(filename = paste0(plotDir, "MeanAnomRich.pdf"), plot = MeanAnomRich, width = 200, height = 150, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("MeanAnomRich.jpeg", device ="jpeg", path = plotDir, width=20, height=15, units="cm", dpi = 350)



##%######################################################%##
#                                                          #
####   Figures, models without interaction with Order   ####
#                                                          #
##%######################################################%##

# load the model outputs if not in environment already
# load(paste0(outDir, "MeanAnomalyModelAbund_noOrder.rdata"))
# load(paste0(outDir, "MeanAnomalyModelRich_noOrder.rdata"))

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

# create table of values for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund2$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund2$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund2$data$LUI)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund2$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund2$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund2$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund2$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund2$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$LUI=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$LUI=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$LUI=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$LUI=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$LUI=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$LUI=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd$LUI <- factor(nd$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


# plot
p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1,1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100), limits = c(-100, 100)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(), 
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("a")


## now the species richness plot 
nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich2$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich2$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich2$data$LUI)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$LUI=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelRich2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich2$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelRich2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich2$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelRich2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich2$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelRich2$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich2$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich2$model,data = nd2, nIters = 10000)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$LUI=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$LUI=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$LUI=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$LUI=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$LUI=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$LUI=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$LUI=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$LUI=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
nd2$LUI <- factor(nd2$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot 
p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_x_continuous(breaks = c(0,0.5, 1,1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125), limits = c(-100, 125)) +
  theme_bw() + 
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  # xlim(c(-1, 5)) +
  # ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("b")


# combine plots
cowplot::plot_grid(p1, p2)

ggsave(filename = paste0(plotDir, "Figure_MeanAnom_Abun_Rich_noOrder.pdf"), plot = last_plot(), width = 183, height = 100, units = "mm", dpi = 300)

