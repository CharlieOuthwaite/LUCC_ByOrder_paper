##%######################################################%##
#                                                          #
####               Land use only models                 ####
#                Tropical vs NonTropical                   #
#                                                          #
##%######################################################%##

# This script takes the processed PREDICTS data, split by latitudinal realm 
# and runs models of the impact of land use and Order on insect biodiversity
# in tropical vs non-tropical realms. 

# clear working environment
rm(list = ls())

# set directories 
inDir<- "1_CheckPrepareData/"
outDir <- "6_TropicalModels/"
plotDir <- "6_TropicalModels/Plots/"
predsDir <- "7_Predictions/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(plotDir)) dir.create(plotDir)
if(!dir.exists(predsDir)) dir.create(predsDir)

# sink(paste0(outDir,"log_LUI_Models_Trop.txt"))
# 
# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
packages_model <- c("StatisticalModels", "predictsFunctions", "lme4", "dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "yarg", "ggplot2", "cowplot", "sjPlot", "gt", "broom.mixed", "MASS")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")

#### 1. Organise data ####

# read in the PREDICTS Site data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds"))

# keep only Coleoptera, Hymenoptera, and Lepidoptera
sites <- filter(sites, Order %in% c('Coleoptera', 'Hymenoptera', 'Lepidoptera'))

# abundance data subset
table(sites[!is.na(sites$LogAbund), 'Realm'])
# NonTropical    Tropical 
#        4674        2270 

# species richness data subset
table(sites[!is.na(sites$Species_richness), 'Realm'])
# NonTropical    Tropical 
#        4887        2405 

# split by land use classes
table(sites$LUI, sites$Realm)
#                      NonTropical Tropical
# Agriculture_High            1478      613
# Agriculture_Low             1014      418
# Primary vegetation          1113      725
# Secondary vegetation        1282      649

# create separate data sets for Tropical and NonTropical sites
trop <- sites[sites$Realm == "Tropical", ]
nontrop <- sites[sites$Realm == "NonTropical", ]

#### 2a. Species richness models, Nontropical ####

# remove NAs in the specified columns
model_data_sr_nontrop <- na.omit(nontrop[,c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
# exclude Neuroptera and Thysanoptera
model_data_sr_nontrop$LUI <- factor(model_data_sr_nontrop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr_nontrop$Order <- factor(model_data_sr_nontrop$Order, levels = c("Coleoptera","Hymenoptera","Lepidoptera"))

# relevel
model_data_sr_nontrop$LUI <- relevel(model_data_sr_nontrop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_sr_nontrop$SS)) # 153
length(unique(model_data_sr_nontrop$SSBS)) # 4040

# look at the spread of land use/use intensity categories
print(table(model_data_sr_nontrop$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               1113                 1282                 1014                 1478 

# effect of land use (the null (intercept-only) model)
sm0_nontrop <-GLMER (modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                     fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of LUI (LUI)
sm3_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                     fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of order only
sm0.2_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                       fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and LUI
sm3.2_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                       fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effect of order and LUI
sm3.3_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                       fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)



#### 2b. Species richness models, Tropical ####

# remove NAs in the specified columns
model_data_sr_trop <- na.omit(trop[,c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_sr_trop$LUI <- factor(model_data_sr_trop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr_trop$Order <- factor(model_data_sr_trop$Order, levels = c("Coleoptera","Hymenoptera","Lepidoptera"))

# relevel
model_data_sr_trop$LUI <- relevel(model_data_sr_trop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_sr_trop$SS)) # 91
length(unique(model_data_sr_trop$SSBS)) # 1683

# look at the spread of land use/use intensity categories
print(table(model_data_sr_trop$LUI)) 

# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#                725                  649                  418                  613 

# Run species richness models using GLMER function from StatisticalModels
# effect of land use (the null (intercept-only) model)
sm0_trop <-GLMER (modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                  fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of LUI (LUI)
sm3_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                  fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of order only
sm0.2_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                    fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and LUI
sm3.2_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                    fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effect of order and LUI
sm3.3_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                    fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)


# take a look at the AICs
AIC_sm_realm<-print(AIC(sm0_trop$model,sm3_trop$model,sm0.2_trop$model,sm3.2_trop$model,sm3.3_trop$model,
                        sm0_nontrop$model,sm3_nontrop$model,sm0.2_nontrop$model,sm3.2_nontrop$model,sm3.3_nontrop$model))

# Warning message:
#   In AIC.default(sm0_trop$model, sm3_trop$model, sm0.2_trop$model,  :
#                    models are not all fitted to the same number of observations

write.csv(AIC_sm_realm,paste0(outDir,"AIC_sm_realm.csv"), row.names = TRUE)

#### 3a. Abundance models, Nontropical ####

# remove NAs in the specified columns
model_data_ab_nontrop <- na.omit(nontrop[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_ab_nontrop$LUI <- factor(model_data_ab_nontrop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab_nontrop$Order <- factor(model_data_ab_nontrop$Order, levels = c("Coleoptera","Hymenoptera","Lepidoptera"))

# relevel
model_data_ab_nontrop$LUI <- relevel(model_data_ab_nontrop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_ab_nontrop$SS)) # 145
length(unique(model_data_ab_nontrop$SSBS)) # 3867

# look at the spread of land use/use intensity categories
print(table(model_data_ab_nontrop$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               1072                 1188                 1004                 1410 

# Run abundance models using 'GLMER' function from StatisticalModels

am0_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                     fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                     fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am0.2_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                       fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.2_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                       fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.3_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                       fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

#### 3b. Abundance models, Tropical ####

# remove NAs in the specified columns
model_data_ab_trop <- na.omit(trop[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_ab_trop$LUI <- factor(model_data_ab_trop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab_trop$Order <- factor(model_data_ab_trop$Order, levels = c("Coleoptera","Hymenoptera","Lepidoptera"))

# relevel
model_data_ab_trop$LUI <- relevel(model_data_ab_trop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_ab_trop$SS)) # 81
length(unique(model_data_ab_trop$SSBS)) # 1548

# look at the spread of land use/use intensity categories
print(table(model_data_ab_trop$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#                674                  606                  411                  579 

# Run abundance models using 'GLMER' function from StatisticalModels

am0_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                  fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                  fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am0.2_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                    fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.2_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                    fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.3_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                    fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

# take a look at the AICs
AIC_am_realm<-print(AIC(am0_trop$model,am3_trop$model,am0.2_trop$model,am3.2_trop$model,am3.3_trop$model,
                        am0_nontrop$model,am3_nontrop$model,am0.2_nontrop$model,am3.2_nontrop$model,am3.3_nontrop$model))

# Warning message:
#   In AIC.default(am0_trop$model, am3_trop$model, am0.2_trop$model,  :
#                    models are not all fitted to the same number of observations

# save as .csv
write.csv(AIC_am_realm, paste0(outDir, "AIC_am_realm.csv"), row.names = TRUE)

# save model data
saveRDS(object = sm0_trop ,file = paste0(outDir,"sm0_trop.rds"))
saveRDS(object = sm3_trop ,file = paste0(outDir,"sm3_trop.rds"))
saveRDS(object = sm0.2_trop ,file = paste0(outDir,"sm0.2_trop.rds"))
saveRDS(object = sm3.2_trop ,file = paste0(outDir,"sm3.2_trop.rds"))
saveRDS(object = sm3.3_trop ,file = paste0(outDir,"sm3.3_trop.rds"))
saveRDS(object = sm0_nontrop ,file = paste0(outDir,"sm0_nontrop.rds"))
saveRDS(object = sm3_nontrop ,file = paste0(outDir,"sm3_nontrop.rds"))
saveRDS(object = sm0.2_nontrop ,file = paste0(outDir,"sm0.2_nontrop.rds"))
saveRDS(object = sm3.2_nontrop ,file = paste0(outDir,"sm3.2_nontrop.rds"))
saveRDS(object = sm3.3_nontrop ,file = paste0(outDir,"sm3.3_nontrop.rds"))
saveRDS(object = am0_trop ,file = paste0(outDir,"am0_trop.rds"))
saveRDS(object = am3_trop ,file = paste0(outDir,"am3_trop.rds"))
saveRDS(object = am0.2_trop ,file = paste0(outDir,"am0.2_trop.rds"))
saveRDS(object = am3.2_trop ,file = paste0(outDir,"am3.2_trop.rds"))
saveRDS(object = am3.3_trop ,file = paste0(outDir,"am3.3_trop.rds"))
saveRDS(object = am0_nontrop ,file = paste0(outDir,"am0_nontrop.rds"))
saveRDS(object = am3_nontrop ,file = paste0(outDir,"am3_nontrop.rds"))
saveRDS(object = am0.2_nontrop ,file = paste0(outDir,"am0.2_nontrop.rds"))
saveRDS(object = am3.2_nontrop ,file = paste0(outDir,"am3.2_nontrop.rds"))
saveRDS(object = am3.3_nontrop ,file = paste0(outDir,"am3.3_nontrop.rds"))
saveRDS(object = model_data_sr_trop ,file = paste0(outDir,"model_data_sr_trop.rds"))
saveRDS(object = model_data_sr_nontrop ,file = paste0(outDir,"model_data_sr_nontrop.rds"))
saveRDS(object = model_data_ab_trop ,file = paste0(outDir,"model_data_ab_trop.rds"))
saveRDS(object = model_data_ab_nontrop ,file = paste0(outDir,"model_data_ab_nontrop.rds"))

##%######################################################%##
#                                                          #
####           Richness and abundance plots             ####
#                 Tropical vs Nontropical                  #
#                                                          #
##%######################################################%##

# read in model data
sm0_trop <- readRDS(file = paste0(outDir,"sm0_trop.rds"))
sm3_trop <- readRDS(file = paste0(outDir,"sm3_trop.rds"))
sm0.2_trop <- readRDS(file = paste0(outDir,"sm0.2_trop.rds"))
sm3.2_trop <- readRDS(file = paste0(outDir,"sm3.2_trop.rds"))
sm3.3_trop <- readRDS(file = paste0(outDir,"sm3.3_trop.rds"))
sm0_nontrop <- readRDS(file = paste0(outDir,"sm0_nontrop.rds"))
sm3_nontrop <- readRDS(file = paste0(outDir,"sm3_nontrop.rds"))
sm0.2_nontrop <- readRDS(file = paste0(outDir,"sm0.2_nontrop.rds"))
sm3.2_nontrop <- readRDS(file = paste0(outDir,"sm3.2_nontrop.rds"))
sm3.3_nontrop <- readRDS(file = paste0(outDir,"sm3.3_nontrop.rds"))
am0_trop <- readRDS(file = paste0(outDir,"am0_trop.rds"))
am3_trop <- readRDS(file = paste0(outDir,"am3_trop.rds"))
am0.2_trop <- readRDS(file = paste0(outDir,"am0.2_trop.rds"))
am3.2_trop <- readRDS(file = paste0(outDir,"am3.2_trop.rds"))
am3.3_trop <- readRDS(file = paste0(outDir,"am3.3_trop.rds"))
am0_nontrop <- readRDS(file = paste0(outDir,"am0_nontrop.rds"))
am3_nontrop <- readRDS(file = paste0(outDir,"am3_nontrop.rds"))
am0.2_nontrop <- readRDS(file = paste0(outDir,"am0.2_nontrop.rds"))
am3.2_nontrop <- readRDS(file = paste0(outDir,"am3.2_nontrop.rds"))
am3.3_nontrop <- readRDS(file = paste0(outDir,"am3.3_nontrop.rds"))
model_data_sr_trop <- readRDS(file = paste0(outDir,"model_data_sr_trop.rds"))
model_data_ab_trop <- readRDS(file = paste0(outDir,"model_data_ab_trop.rds"))
model_data_sr_nontrop <- readRDS(file = paste0(outDir,"model_data_sr_nontrop.rds"))
model_data_ab_nontrop <- readRDS(file = paste0(outDir,"model_data_ab_nontrop.rds"))

#### 4. Model selection tables ####


# table of AICs
selection_table_trop <- data.frame("Realm" = c(rep("Tropical", 10)),
                                   "Response" = c(rep("Species richness", 5),
                                                  rep("Total abundance", 5)),
                                   "Model" = c("Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ Order + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                                   "AIC" = c(AIC(sm0_trop$model), AIC(sm3_trop$model), AIC(sm0.2_trop$model), AIC(sm3.2_trop$model), AIC(sm3.3_trop$model),
                                             AIC(am0_trop$model), AIC(am3_trop$model), AIC(am0.2_trop$model), AIC(am3.2_trop$model),AIC(am3.3_trop$model))) %>%
  group_by(Response) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

selection_table_nontrop <- data.frame("Realm" = c(rep("NonTropical", 10)),
                                      "Response" = c(rep("Species richness", 5),
                                                     rep("Total abundance", 5)),
                                      "Model" = c("Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ Order + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                                      "AIC" = c(AIC(sm0_nontrop$model), AIC(sm3_nontrop$model), AIC(sm0.2_nontrop$model), AIC(sm3.2_nontrop$model), AIC(sm3.3_nontrop$model),  
                                                AIC(am0_nontrop$model), AIC(am3_nontrop$model), AIC(am0.2_nontrop$model), AIC(am3.2_nontrop$model),AIC(am3.3_nontrop$model))) %>%
  group_by(Response) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# save the tables
gtsave(selection_table_trop,paste0(outDir,"SimpleLUIModel_Selection_trop.html"))
gtsave(selection_table_trop,"SimpleLUIModel_Selection_trop.png", path = outDir)

gtsave(selection_table_nontrop,paste0(outDir,"SimpleLUIModel_Selection_nontrop.html"))
gtsave(selection_table_nontrop,"SimpleLUIModel_Selection_nontrop.png", path = outDir)

####  5a. Species Richness Plot, Nontropical ####

richness_metric_nontrop <- predict_effects(iterations = 1000,
                                           model = sm3.3_nontrop$model,
                                           model_data = model_data_sr_nontrop,
                                           response_variable = "Species_richness",
                                           fixed_number = 2,
                                           fixed_column = c("Order", "LUI"),
                                           factor_number_1 = 3,
                                           factor_number_2 = 4,
                                           neg_binom = FALSE)
# predictions
result.sr.nontrop <- fin_conf
result.sr.nontrop <- dplyr::select(result.sr.nontrop,-c(Species_richness))
richness_metric_nontrop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_sr_nontrop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric_nontrop)

####  5b. Species Richness Plot, Tropical ####

richness_metric_trop <- predict_effects(iterations = 1000,
                                        model = sm3.3_trop$model,
                                        model_data = model_data_sr_trop,
                                        response_variable = "Species_richness",
                                        fixed_number = 2,
                                        fixed_column = c("Order", "LUI"),
                                        factor_number_1 = 3,
                                        factor_number_2 = 4,
                                        neg_binom = FALSE)
result.sr.trop <- fin_conf
result.sr.trop <- dplyr::select(result.sr.trop,-c(Species_richness))
richness_metric_trop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_sr_trop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric_trop)

# plot both realms together
richness_nontrop <- richness_metric_nontrop + 
  labs(y ="Species richness diff. (%)", x = "Order") +
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 150)) +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "none") 

richness_trop <- richness_metric_trop +
  labs(y ="Species richness diff. (%)", x = "Order") +
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 150)) +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "none")

legend <- get_legend(
  richness_nontrop +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.box = "horizontal", 
          legend.text = element_text(size = 7), 
          legend.title = element_blank())
)

# put together
richness_realms <- cowplot::plot_grid(NULL,NULL,richness_nontrop,richness_trop, ncol = 2, rel_heights = c(0.1,1),labels=c('','','a - Non-tropical','b - Tropical'),
                                      label_size = 10,
                                      label_x = 0, label_y = 1,
                                      hjust = -0.2, vjust = -0.1)
# add legend
richness_realms_LUI <- cowplot::plot_grid(richness_realms, legend, ncol = 1, rel_heights = c(1,0.1))

# save plot (pdf)
ggsave(filename = paste0(plotDir, "SimpleLUI_Rich_Realms.pdf"), plot = richness_realms_LUI, width = 200, height = 100, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("SimpleLUI_Rich_Realms.jpeg", device ="jpeg", path = plotDir, width=200, height=100, units="mm", dpi = 350)

####  6a. Abundance Plot, Nontropical ####

abundance_metric_nontrop <- predict_effects(iterations = 1000,
                                            model = am3.3_nontrop$model,
                                            model_data = model_data_ab_nontrop,
                                            response_variable = "LogAbund",
                                            fixed_number = 2,
                                            fixed_column = c("Order", "LUI"),
                                            factor_number_1 = 3,
                                            factor_number_2 = 4,
                                            neg_binom = FALSE)
result.ab.nontrop <- fin_conf
result.ab.nontrop <- dplyr::select(result.ab.nontrop,-c(LogAbund))
abundance_metric_nontrop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab_nontrop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric_nontrop)

####  6b. Abundance Plot, Tropical ####

abundance_metric_trop <- predict_effects(iterations = 1000,
                                         model = am3.3_trop$model,
                                         model_data = model_data_ab_trop,
                                         response_variable = "LogAbund",
                                         fixed_number = 2,
                                         fixed_column = c("Order", "LUI"),
                                         factor_number_1 = 3,
                                         factor_number_2 = 4,
                                         neg_binom = FALSE)
result.ab.trop <- fin_conf
result.ab.trop <- dplyr::select(result.ab.trop,-c(LogAbund))
abundance_metric_trop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab_trop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric_trop)

# plot both realms together
abundance_nontrop <- abundance_metric_nontrop + 
  labs(y ="Total abundance diff. (%)", x = "Order") +
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 150)) +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "none")

abundance_trop <- abundance_metric_trop +
  labs(y ="Total abundance diff. (%)", x = "Order") +
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 150)) +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "none")

legend <- get_legend(
  richness_nontrop +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.box = "horizontal", 
          legend.text = element_text(size = 7), 
          legend.title = element_blank())
)

# put together
abundance_realms <- cowplot::plot_grid(NULL,NULL,abundance_nontrop,abundance_trop, ncol = 2, rel_heights = c(0.1,1),labels=c('','','a - Non-tropical','b - Tropical'),
                                      label_size = 10,
                                      label_x = 0, label_y = 1,
                                      hjust = -0.2, vjust = -0.1)
# add legend
abundance_realms_LUI <- cowplot::plot_grid(abundance_realms, legend, ncol = 1, rel_heights = c(1,0.1))

# save plot (pdf)
ggsave(filename = paste0(plotDir, "SimpleLUI_Abund_Realms.pdf"), plot = abundance_realms_LUI, width = 200, height = 100, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("SimpleLUI_Abund_Realms.jpeg", device ="jpeg", path = plotDir, width=200, height=100, units="mm", dpi = 350)

####  7. Plot all together  ####

# plot species richness and abundance responses from both realms together #
simplemods_realms<-cowplot::plot_grid(NULL,NULL,richness_nontrop,richness_trop,NULL,NULL,abundance_nontrop,abundance_trop, ncol=2, rel_heights = c(0.1,1,0.1,1),labels=c('','','a - Non-tropical','b - Tropical','','','c - Non-tropical','d - Tropical'),
                                      label_size = 10,
                                      label_x = 0, label_y = 1,
                                      hjust = -0.2, vjust = -0.1)
# add legend
simplemods_realms <- cowplot::plot_grid(simplemods_realms, legend, ncol = 1, rel_heights = c(1,0.1))

# save plot (pdf)
ggsave(filename = paste0(plotDir, "SimpleLUI_Realms.pdf"), plot = simplemods_realms, width = 200, height = 200, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("SimpleLUI_Realms.jpeg", device ="jpeg", path = plotDir, width=200, height=200, units="mm", dpi = 350)

#### 8. Table of predicted values ####

# combine results into a table for saving
all_res <- rbind(result.ab.nontrop, result.ab.trop, result.sr.nontrop, result.sr.trop)
all_res$measure <- c(rep("ab",24 ), rep("sr", 24))
all_res$Realm <- c(rep("nontrop",12), rep("trop", 12),rep("nontrop",12), rep("trop",12))

# save .csv
write.csv(all_res, file = paste0(predsDir,"LUI_predictions_realms.csv"))

# save as gt table
LUI_predictions_realms <- all_res %>% gt()
gtsave(LUI_predictions_realms,"LUI_predictions_realms.png", path = predsDir)



# R-squared values

R2GLMER(am3.3_trop$model)
# $conditional
# [1] 0.7149983
# 
# $marginal
# [1] 0.06486359

R2GLMER(am3.3_nontrop$model)
# $conditional
# [1] 0.8209407
# 
# $marginal
# [1] 0.09744095

R2GLMER(sm3.3_trop$model)
# $conditional
# [1] 0.7040361
# 
# $marginal
# [1] 0.0782471

R2GLMER(sm3.3_nontrop$model)
# $conditional
# [1] 0.7436253
# 
# $marginal
# [1] 0.07648085


# t.end <- Sys.time()
# 
# print(round(t.end - t.start,0))
# 
# sink()

