##%######################################################%##
#                                                          #
####    Run models for climate land use interactions    ####
#                Tropical vs NonTropical                   #
#                                                          #
##%######################################################%##

# This script takes the processed PREDICTS data, split by latitudinal realm 
# and models are run to assess the response of insect 
# biodiversity to the land use intensity metric and to the climate
# anomalies and their interaction.

# clear environment
rm(list = ls())

# set directories 
predictsDir<- "5_RunLUIClimateModels/"
outDir <- "6_TropicalModels/"
if(!dir.exists(outDir)) dir.create(outDir)

# sink(paste0(outDir,"log_LUI_ClimateModels.txt"))
# 
# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
packages_model <- c("StatisticalModels", "predictsFunctions", "lme4", "dplyr","devtools")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

# source additional functions
source("0_Functions.R")

#### 1. Organise data ####

# read in the predicts data from the LUI Climate Models
predictsSites <- readRDS(paste0(predictsDir,"PREDICTSSitesClimate_Data.rds"))

# keep only Coleoptera, Hymenoptera, and Lepidoptera
predictsSites <- filter(predictsSites, Order %in% c('Coleoptera', 'Hymenoptera', 'Lepidoptera'))

# split into two data frames based on latitudinal realm

trop <- predictsSites[predictsSites$Realm == "Tropical", ]
nontrop <- predictsSites[predictsSites$Realm == "NonTropical", ]

# take a look at possible correlations between variables
# Tropical
cor(trop$avg_temp, trop$TmeanAnomaly)

# 0.01349527

cor(trop$avg_temp, trop$StdTmeanAnomaly)

# 0.2865458

cor(trop$TmeanAnomaly, trop$StdTmeanAnomaly)

# 0.2865458

# NonTropical
cor(nontrop$avg_temp, nontrop$TmeanAnomaly)

# -0.112796

cor(nontrop$avg_temp, nontrop$StdTmeanAnomaly)

# -0.6327842

cor(nontrop$TmeanAnomaly, nontrop$StdTmeanAnomaly)

# 0.6392614

# save the datasets
saveRDS(object = trop,file = paste0(outDir,"trop.rds"))
saveRDS(object = nontrop,file = paste0(outDir,"nontrop.rds"))


#### 2. Create models of climate land use interactions ####

# ia. Abundance, Mean Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund_nontrop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                       fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                                       randomStruct = "(1|SS)+(1|SSB)",
                                       saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund_nontrop$model)

# save the model output
save(MeanAnomalyModelAbund_nontrop, file = paste0(outDir, "MeanAnomalyModelAbund_nontrop.rdata"))

# iia. Richness, Mean Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich_nontrop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                      fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                                      randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                      saveVars = c("SSBS"))

summary(MeanAnomalyModelRich_nontrop$model)

save(MeanAnomalyModelRich_nontrop, file = paste0(outDir, "MeanAnomalyModelRich_nontrop.rdata"))

# iiia. Abundance, Max Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrops$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund_nontrop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                      fixedStruct = "LUI * StdTmaxAnomalyRS * Order",
                                      randomStruct = "(1|SS)+(1|SSB)",
                                      saveVars = c("SSBS"))

summary(MaxAnomalyModelAbund_nontrop$model)

save(MaxAnomalyModelAbund_nontrop, file = paste0(outDir, "MaxAnomalyModelAbund_nontrop.rdata"))

# iva. Richness, Max Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$StdTmeanAnomalyRS),]

MaxAnomalyModelRich_nontrop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                     fixedStruct = "LUI * StdTmaxAnomalyRS * Order",
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                     saveVars = c("SSBS"))

# save model output
save(MaxAnomalyModelRich_nontrop, file = paste0(outDir, "MaxAnomalyModelRich_nontrop.rdata"))

# ib. Abundance, mean anomaly, Tropical

model_data <- trop[!is.na(trop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund_trop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                    fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                                    randomStruct = "(1|SS)+(1|SSB)",
                                    saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund_trop$model)

# save the model output
save(MeanAnomalyModelAbund_trop, file = paste0(outDir, "MeanAnomalyModelAbund_trop.rdata"))

# iib. Richness, mean anomaly, Tropical

model_data <- trop[!is.na(trop$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich_trop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                   fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                   saveVars = c("SSBS"))

summary(MeanAnomalyModelRich_trop$model)

save(MeanAnomalyModelRich_trop, file = paste0(outDir, "MeanAnomalyModelRich_trop.rdata"))

# iiib. Abundance, Max Anomaly, Tropical

model_data <- trop[!is.na(trop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund_trop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                   fixedStruct = "LUI * StdTmaxAnomalyRS * Order",
                                   randomStruct = "(1|SS)+(1|SSB)",
                                   saveVars = c("SSBS"))

summary(MaxAnomalyModelAbund_trop$model)

save(MaxAnomalyModelAbund_trop, file = paste0(outDir, "MaxAnomalyModelAbund_trop.rdata"))

# ivb. Richness, Max Anomaly, Tropical

model_data <- trop[!is.na(trop$StdTmeanAnomalyRS),]

MaxAnomalyModelRich_trop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                  fixedStruct = "LUI * StdTmaxAnomalyRS * Order",
                                  randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                  saveVars = c("SSBS"))

summary(MaxAnomalyModelRich_trop$model)

# save model output
save(MaxAnomalyModelRich_trop, file = paste0(outDir, "MaxAnomalyModelRich_trop.rdata"))


#### 3. Model output tables ####

# nontropical

tab_model(MeanAnomalyModelAbund_nontrop$model, transform = NULL, file = paste0(tabDir,"AbunMeanAnom_output_nontrop_table.html"))
summary(MeanAnomalyModelAbund_nontrop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_nontrop$model) # check the R2 values 
# $conditional
# [1] 0.5761286
# 
# $marginal
# [1] 0.1504236

tab_model(MeanAnomalyModelRich_nontrop$model, transform = NULL, file = paste0(tabDir,"RichMeanAnom_output_nontrop_table.html"))
summary(MeanAnomalyModelRich_nontrop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_nontrop$model) # check the R2 values
# $conditional
# [1] 0.7402805
# 
# $marginal
# [1] 0.1029864

tab_model(MaxAnomalyModelAbund_nontrop$model, transform = NULL, file = paste0(tabDir,"AbunMaxAnom_output_nontrop_table.html"))
summary(MaxAnomalyModelAbund_nontrop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_nontrop$model) # check the R2 values 

tab_model(MaxAnomalyModelRich_nontrop$model, transform = NULL, file = paste0(tabDir,"RichMaxAnom_output_nontrop_table.html"))
summary(MaxAnomalyModelRich_nontrop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_nontrop$model) # check the R2 values 

# tropical

tab_model(MeanAnomalyModelAbund_trop$model, transform = NULL, file = paste0(tabDir,"AbunMeanAnom_output_trop_table.html"))
summary(MeanAnomalyModelAbund_trop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_trop$model) # check the R2 values 
# $conditional
# [1] 0.4410931
# 
# $marginal
# [1] 0.1449011

tab_model(MeanAnomalyModelRich_trop$model, transform = NULL, file = paste0(tabDir, "RichMeanAnom_output_trop_table.html"))
summary(MeanAnomalyModelRich_trop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_trop$model) # check the R2 values
# $conditional
# [1] 0.7506744
# 
# $marginal
# [1] 0.1304144


tab_model(MaxAnomalyModelAbund_trop$model, transform = NULL, file = paste0(tabDir,"AbunMaxAnom_output_trop_table.html"))
summary(MaxAnomalyModelAbund_trop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_trop$model) # check the R2 values 

tab_model(MaxAnomalyModelRich_trop$model, transform = NULL, file = paste0(tabDir,"RichMaxAnom_output_trop_table.html"))
summary(MaxAnomalyModelRich_trop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_trop$model) # check the R2 values 


