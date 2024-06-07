##%######################################################%##
#                                                          #
####                  Bayesian models                   ####
#                                                          #
##%######################################################%##

# Here I run the same model formulations but in a Bayesian framework


# ensure working directory is clear
rm(list = ls())



#### 1. Land use only models ####

# set up directories
inDir <- "1_CheckPrepareData/"
outDir <- "2_RunSimpleLUIModel/"
predsDir <- "7_Predictions/"
dataDir <- "Data/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(predsDir)) dir.create(predsDir)


# load datasets

# load richness model data, model_data_sr
saveRDS(object = model_data_sr ,file = paste0(outDir,"model_data_sr.rds"))
# load abundance model data, model_data_ab
saveRDS(object = model_data_ab ,file = paste0(outDir,"model_data_ab.rds"))


#### 2. Land use and STA models ####