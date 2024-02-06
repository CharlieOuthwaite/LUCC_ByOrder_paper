##%######################################################%##
#                                                          #
####                    MCMC models                     ####
#                                                          #
##%######################################################%##

# In this script the models are rerun using a Bayesian approach. 

# load packages
library(brms)
library(sjPlot)

# directories
datadir1 <- "2_RunSimpleLUIModel/"
datadir2 <- "5_RunLUIClimateModels/"
outdir <- "10_Additional_Tests/MCMCModels/"
if(!dir.exists(outdir)) dir.create(outdir)

# load in the datasets
sr_data <- readRDS(paste0(datadir1, "model_data_sr.rds")) # LU only model
ab_data <- readRDS(paste0(datadir1, "model_data_ab.rds")) # LU only model

pred_data <- readRDS(paste0(datadir2, "PREDICTSSitesClimate_Data.rds")) # LU STA model


############################################################
#                                                          #
####              1. Land use only models               ####
#                                                          #
############################################################

# 
# ## richness
# rich_LU <- brm(Species_richness ~ LUI + (1|SS)+(1|SSB)+(1|SSBS), 
#                data = sr_data, 
#                family = poisson("log"), 
#                chains = 4, 
#                iter = 5000, 
#                warmup = 1000, 
#                thin = 1)
# 
# # take a look
# rich_LU
# # check plots
# plot(rich_LU)
# # check model fit
# pp_check(rich_LU)
# # save output
# save(rich_LU, file = paste0(outdir, "Richness_LU_bayesmod.Rdata"))
# 
# 
# ## abundance
# abun_LU <- brm(LogAbund ~ LUI + (1|SS)+(1|SSB), 
#                data = ab_data, 
#                family = gaussian(), # also ran skew_normal()
#                chains = 4, 
#                iter = 5000, 
#                warmup = 1000, 
#                thin = 1)
# 
# # take a look
# abun_LU
# # check plots
# plot(abun_LU)
# # check model fit
# pp_check(abun_LU)
# # save output
# save(abun_LU, file = paste0(outdir, "Abundance_LU_bayesmod.Rdata"))


############################################################
#                                                          #
####             2. Land use and Order                  ####
#                                                          #
############################################################

#https://ourcodingclub.github.io/tutorials/brms/

## richness
rich_LUOr <- brm(Species_richness ~ Order * LUI + (1|SS)+(1|SSB)+(1|SSBS), 
               data = sr_data, 
               family = poisson("log"), 
               cores = 4)

# take a look
rich_LUOr
# check plots
plot(rich_LUOr)
# check model fit
pp_check(rich_LUOr)
# save output
save(rich_LUOr, file = paste0(outdir, "Richness_LUOrder_bayesmod.Rdata"))



## abundance
abun_LUOr <- brm(LogAbund ~ Order * LUI + (1|SS)+(1|SSB), 
               data = ab_data, 
               family = gaussian(), 
               cores = 4)

# take a look
abun_LUOr
# check plots
plot(abun_LUOr)
# check model fit
pp_check(abun_LUOr)
# save output
save(abun_LUOr, file = paste0(outdir, "Abundance_LUOrder_bayesmod.Rdata"))


# 
# # 3. Land use and climate ------------------------------------------------------
# 
# 
# ## richness
# rich_LUSTA <- brm(Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS)+(1|SSB)+(1|SSBS), 
#                  data = pred_data, 
#                  family = poisson("log"), 
#                  chains = 4, 
#                  iter = 5000, 
#                  warmup = 1000, 
#                  thin = 1)
# 
# # take a look
# rich_LUSTA
# # check plots
# plot(rich_LUSTA)
# # check model fit
# pp_check(rich_LUSTA)
# # save output
# save(rich_LUSTA, file = paste0(outdir, "Richness_LUSTA_bayesmod.Rdata"))
# 
# 
# 
# # subset to abundance data
# pred_abun <- pred_data[!is.na(pred_data$LogAbund), ]
# pred_abun <- droplevels(pred_abun)
# 
# ## abundance
# abun_LUSTA <- brm(LogAbund ~ LUI * StdTmeanAnomalyRS + (1|SS)+(1|SSB), 
#                  data = pred_abun, 
#                  family = gaussian(), 
#                  chains = 4, 
#                  iter = 5000, 
#                  warmup = 1000, 
#                  thin = 1)
# 
# # take a look
# abun_LUSTA
# # check plots
# plot(abun_LUSTA)
# # check model fit
# pp_check(abun_LUSTA)
# # save output
# save(abun_LUSTA, file = paste0(outdir, "Abundance_LUSTA_bayesmod.Rdata"))
# 


# 4. Land use, climate and order------------------------------------------------

## richness
rich_LUSTAOr <- brm(Species_richness ~ LUI * StdTmeanAnomalyRS * Order + (1|SS)+(1|SSB)+(1|SSBS), 
                  data = pred_data, 
                  family = poisson("log"), 
                  cores = 4)

# take a look
rich_LUSTAOr
# check plots
plot(rich_LUSTAOr)
# check model fit
pp_check(rich_LUSTAOr)
# save output
save(rich_LUSTAOr, file = paste0(outdir, "Richness_LUSTAOr_bayesmod.Rdata"))



# subset to the abundance data
pred_abun <- pred_data[!is.na(pred_data$LogAbund), ]

## abundance
abun_LUSTAOr <- brm(LogAbund ~ LUI * StdTmeanAnomalyRS * Order + (1|SS)+(1|SSB), 
                  data = pred_abun, 
                  family = gaussian(), 
                  cores = 4)

# take a look
abun_LUSTAOr
# check plots
plot(abun_LUSTAOr)
# check model fit
pp_check(abun_LUSTAOr)
# save output
save(abun_LUSTAOr, file = paste0(outdir, "Abundance_LUSTAOr_bayesmod.Rdata"))




############################################################
#                                                          #
#              Save model outputs as tables                #
#                                                          #
############################################################


# LU only models, richness
tab_model(rich_LUOr, transform = NULL, file = paste0(outdir,"Bayes_rich_LU_output_table.html"), show.icc = F)
# LU only models, abundance
tab_model(abun_LUOr, transform = NULL, file = paste0(outdir,"Bayes_abun_LU_output_table.html"), show.icc = F)
# LU and STA models, richness
tab_model(rich_LUSTAOr, transform = NULL, file = paste0(outdir,"Bayes_rich_LUSTA_output_table.html"), show.icc = F)
# LU and STA models, abundance
tab_model(abun_LUSTAOr, transform = NULL, file = paste0(outdir,"Bayes_abun_LUSTA_output_table.html"), show.icc = F)




############################################################
#                                                          #
#                     plot the results                     #
#                                                          #
############################################################






############################################################
#                                                          #
#              Run separate models per Order               #
#                                                          #
############################################################




#########  Land use and climate  ###########

# load in the data
predictsSites <- readRDS(file = paste0(datadir,"PREDICTSSitesClimate_Data.rds"))


#### loop through and run a model for each order ####

orders <- c("Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera")

# order <- orders[5]
for(order in orders){
  
  print(order)
  
  # subset the data to just those sites for the order of interest
  order_data <- predictsSites[predictsSites$Order == order, ]
  
  # remove any rows with NAs in 
  abun_data <- order_data[!is.na(order_data$LogAbund), ]
  
  print("abundance")
  
  ## abundance
  order_mod_ab   <- brm(LogAbund ~ LUI * StdTmeanAnomalyRS + (1|SS)+(1|SSB), 
                     data = abun_data, 
                     family = gaussian(), 
                     cores = 4)
  
  
  # save the model output
  save(order_mod_ab, file = paste0(outdir, "MeanAnomalyModelAbund_", order, ".rdata"))
  
  
  ## richness
  
  print("richness")

  order_mod_sr <- brm(Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS)+(1|SSB)+(1|SSBS), 
                      data = pred_data, 
                      family = poisson("log"), 
                      cores = 4)
  
  save(order_mod_sr, file = paste0(outdir, "MeanAnomalyModelRich_", order, ".rdata"))
  
  
  
}


