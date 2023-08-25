##%######################################################%##
#                                                          #
####                    MCMC models                     ####
#                                                          #
##%######################################################%##

# In this script the models are rerun using a Bayesian approach. 

# load packages
library(brms)

# directories
datadir1 <- "2_RunSimpleLUIModel/"
datadir2 <- "5_RunLUIClimateModels/"

# load in the datasets
sr_data <- readRDS(paste0(datadir1, "model_data_sr.rds")) # LU only model
ab_data <- readRDS(paste0(datadir1, "model_data_ab.rds")) # LU only model

pred_data <- readRDS(paste0(datadir2, "PREDICTSSitesClimate_Data.rds")) # LU STA model


############################################################
#                                                          #
#                   Land use only models                   #
#                                                          #
############################################################

# 1. Land use only -------------------------------------------------------------

## richness
rich_LU <- brm(Species_richness ~ LUI + (1|SS)+(1|SSB)+(1|SSBS), 
               data = sr_data, 
               family = poisson("log"), 
               chains = 4, 
               iter = 5000, 
               warmup = 1000, 
               thin = 1)



## abundance
rich_LU <- brm(LogAbund ~ LUI + (1|SS)+(1|SSB), 
               data = ab_data, 
               family = gaussian(), 
               chains = 4, 
               iter = 5000, 
               warmup = 1000, 
               thin = 1)

# 2. Land use and Order --------------------------------------------------------

## richness
rich_LU <- brm(Species_richness ~ LUI * Order + (1|SS)+(1|SSB)+(1|SSBS), 
               data = sr_data, 
               family = poisson("log"), 
               chains = 4, 
               iter = 5000, 
               warmup = 1000, 
               thin = 1)



## abundance
rich_LU <- brm(LogAbund ~ LUI * Order + (1|SS)+(1|SSB), 
               data = ab_data, 
               family = gaussian(), 
               chains = 4, 
               iter = 5000, 
               warmup = 1000, 
               thin = 1)






