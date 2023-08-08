############################################################
#                                                          #
#                   Land use only models                   #
#                         Global                           #
#                                                          #
############################################################

# This script takes the processed PREDICTS data and runs models of the 
# impact of land use and Order on insect biodiversity. 


# ensure working directory is clear
rm(list = ls())

# set up directories
inDir <- "1_CheckPrepareData/"
outDir <- "2_RunSimpleLUIModel/"
predsDir <- "7_Predictions/"
dataDir <- "Data/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(predsDir)) dir.create(predsDir)

# load libraries
packages_model <- c("StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS","webshot")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")


#### 1. Organise data ####

# read in the Site data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds")) # 8884 rows

## Species Richness Model ##

# remove NAs in the specified columns
model_data_sr <- na.omit(sites[,c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')]) # 8884 rows

# order data
model_data_sr$LUI <- factor(model_data_sr$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr$Order <- factor(model_data_sr$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"))

#save model_data_sr
saveRDS(object = model_data_sr ,file = paste0(outDir,"model_data_sr.rds"))

# summaries
length(unique(model_data_sr$SS)) # 254
length(unique(model_data_sr$SSBS)) # 6014

table(model_data_sr$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
# 2315         721         871        3328        1649

# look at the spread of land use/use intensity categories
print(table(model_data_sr$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               2093                 2392                 1886                 2513  


#### 2. Species Richness models ####

# Run species richness models using GLMER function from StatisticalModels

# null (intercept-only) model
sm0 <-GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
            fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of Land Use
sm1 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of land use and use intensity
sm2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of LUI
sm3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# check the R2 values
R2GLMER(sm3$model)
# $conditional
# [1] 0.6700198
# 
# $marginal
# [1] 0.007716616

# interactive effects of land use and use intensity
sm4 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of order only
sm0.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and land use
sm1.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and land use and use intensity
sm2.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and LUI
sm3.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and interactive effects of land use and use intensity
sm4.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# interactive effect of order and land use
sm1.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# interactive effects of order and additive effects of land use and use intensity
sm2.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effect of order and LUI
sm3.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# check the R2 values
R2GLMER(sm3.3$model) # with order level interaction
# $conditional
# [1] 0.7106989
# 
# $marginal
# [1] 0.05387385

# interactive effects of order, land use, and use intensity
sm4.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# take a look at the AICs
AIC_sm<-print(AIC(sm0$model,sm1$model,sm2$model,sm3$model,sm4$model,
                  sm0.2$model,sm1.2$model,sm2.2$model,sm3.2$model,sm4.2$model,
                  sm1.3$model,sm2.3$model,sm3.3$model,sm4.3$model))

write.csv(AIC_sm, file = paste0(outDir,"AIC_sm.csv"))

# Run richness models using 'glmr' function from lme4
# need these to run allFit()
# need to run allFit() to test the models

g_sm0 <- glmer(Species_richness~1+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm1 <- glmer(Species_richness~LandUse+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm2 <- glmer(Species_richness~LandUse+Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm3 <- glmer(Species_richness~LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm4 <- glmer(Species_richness~LandUse*Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm0.2 <- glmer(Species_richness~Order+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm1.2 <- glmer(Species_richness~Order+LandUse+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm2.2 <- glmer(Species_richness~Order+LandUse+Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm3.2 <- glmer(Species_richness~Order+LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm4.2 <- glmer(Species_richness~Order+LandUse*Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm1.3 <- glmer(Species_richness~Order*LandUse+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm2.3 <- glmer(Species_richness~Order*LandUse+Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm3.3 <- glmer(Species_richness~Order*LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm4.3 <- glmer(Species_richness~Order*LandUse*Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)

# take a look at the AICs
AIC_g_sm<-print(AIC(g_sm0,g_sm1,g_sm2,g_sm3,g_sm4,
                    g_sm0.2,g_sm1.2,g_sm2.2,g_sm3.2,g_sm4.2,
                    g_sm1.3,g_sm2.3,g_sm3.3,g_sm4.3))

write.csv(AIC_g_sm, file = paste0(outDir,"AIC_g_sm.csv"))

# Testing richness models with allFit()

g_sm0.2_all <-allFit(g_sm0.2)
g_sm1.2_all <-allFit(g_sm1.2)
g_sm2.2_all <-allFit(g_sm2.2)
g_sm3.2_all <-allFit(g_sm3.2)
g_sm4.2_all <-allFit(g_sm4.2)
g_sm1.3_all <-allFit(g_sm1.3)
g_sm2.3_all <-allFit(g_sm2.3)
g_sm3.3_all <-allFit(g_sm3.3)
g_sm4.3_all <-allFit(g_sm4.3)

save(sm3.3, file = paste0(outDir, "Richness_landuse_model.rdata"))

#### 3. Abundance Models  #####

# remove NAs in the specified columns
model_data_ab <- na.omit(sites[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_ab$LUI <- factor(model_data_ab$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab$Order <- factor(model_data_ab$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"))

# relevel
model_data_ab$LUI <- relevel(model_data_ab$LUI, ref = "Primary vegetation")

#save model_data_ab
saveRDS(object = model_data_ab ,file = paste0(outDir,"model_data_ab.rds"))

# summaries
length(unique(model_data_ab$SS)) # 236 Studies
length(unique(model_data_ab$SSBS))# 5706

table(model_data_ab$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       2295         701         851        3146        1503 

# look at the spread of land use/use intensity categories
print(table(model_data_ab$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               2001                 2255                 1869                 2371

# Run abundance models using 'GLMER' function from StatisticalModels

am0 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am1 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

# check the R2 values
R2GLMER(am3$model)
# $conditional
# [1] 0.6838617
# 
# $marginal
# [1] 0.007807588

am4 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am0.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am1.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order+LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am2.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order+LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.2<- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
              fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am4.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order+LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am1.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am2.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

# check the R2 values
R2GLMER(am3.3$model)
# $conditional
# [1] 0.7623785
# 
# $marginal
# [1] 0.06004844

am3.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am4.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

# take a look at the AICs
AIC_am <-print(AIC(am0$model,am1$model,am2$model,am3$model,am4$model,
                   am0.2$model,am1.2$model,am2.2$model,am3.2$model,am4.2$model,
                   am1.3$model,am2.3$model,am3.3$model,am4.3$model))

write.csv(AIC_am, file = paste0(outDir,"AIC_am.csv"))

# save 
saveRDS(object = sm0 ,file = paste0(outDir,"sm0.rd"))
saveRDS(object = sm3 ,file = paste0(outDir,"sm3.rds"))
saveRDS(object = sm0.2 ,file = paste0(outDir,"sm0.2.rds"))
saveRDS(object = sm3.2 ,file = paste0(outDir,"sm3.2.rds"))
saveRDS(object = sm3.3 ,file = paste0(outDir,"sm3.3.rds"))
saveRDS(object = am0 ,file = paste0(outDir,"am0.rds"))
saveRDS(object = am3 ,file = paste0(outDir,"am3.rds"))
saveRDS(object = am0.2 ,file = paste0(outDir,"am0.2.rds"))
saveRDS(object = am3.2 ,file = paste0(outDir,"am3.2.rds"))
saveRDS(object = am3.3 ,file = paste0(outDir,"am3.3.rds"))

save(am3.3, file = paste0(outDir, "Abundance_landuse_model.rdata"))


# read in model data
sm0 <- readRDS(file = paste0(outDir,"sm0.rds"))
sm3 <- readRDS(file = paste0(outDir,"sm3.rds"))
sm0.2 <- readRDS(file = paste0(outDir,"sm0.2.rds"))
sm3.2 <- readRDS(file = paste0(outDir,"sm3.2.rds"))
sm3.3 <- readRDS(file = paste0(outDir,"sm3.3.rds"))
am0 <- readRDS(file = paste0(outDir,"am0.rds"))
am3 <- readRDS(file = paste0(outDir,"am3.rds"))
am0.2 <- readRDS(file = paste0(outDir,"am0.2.rds"))
am3.2 <- readRDS(file = paste0(outDir,"am3.2.rds"))
am3.3 <- readRDS(file = paste0(outDir,"am3.3.rds"))
model_data_sr <- readRDS(file = paste0(outDir,"model_data_sr.rds"))
model_data_ab <- readRDS(file = paste0(outDir,"model_data_ab.rds"))


#### 4. Table of AICs ####

# species richness and abundance together
selection_table <- data.frame("Response" = c(rep("Species richness", 5),
                                             rep("Abundance", 5)),
                              "Model" = c("Species richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Abundance ~ 1 + (1|SS) + (1|SSB)",
                                          "Abundance ~ LUI + (1|SS) + (1|SSB)",
                                          "Abundance ~ Order + (1|SS) + (1|SSB)",
                                          "Abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                          "Abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                              "AIC" = c(AIC(sm0$model), AIC(sm3$model), AIC(sm0.2$model), AIC(sm3.2$model), AIC(sm3.3$model),  
                                        AIC(am0$model), AIC(am3$model), AIC(am0.2$model), AIC(am3.2$model),AIC(am3.3$model))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  dplyr::select(Model,AIC,deltaAIC) %>%
  gt(rowname_col = "Model") %>%
  tab_row_group(
    label = "Species Richness",
    rows = starts_with("Species richness")
  ) %>%
  tab_row_group(
    label = "Abundance",
    rows = starts_with("Abundance")
  )%>% 
  cols_align(
    align = "center",
    columns = c(Model, AIC, deltaAIC)
  )%>%
  tab_stubhead(label = "Models")

# save it
gtsave(selection_table, "SimpleLUIModel_Selection.png", path = outDir)

# species richness only
selection_table_sr <- data.frame("Response" = c(rep("Species richness", 5)),
                                 "Model" = c("Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)"),
                                 "AIC" = c(AIC(sm0$model), AIC(sm3$model), AIC(sm0.2$model), AIC(sm3.2$model), AIC(sm3.3$model))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# save it
gtsave(selection_table_sr,"SimpleLUIModel_Selection_Rich.png", path = outDir)

# total abundance only
selection_table_ab <- data.frame("Response" = c(rep("Total abundance", 5)),
                                 "Model" = c("Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ Order + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                                 "AIC" = c(AIC(am0$model), AIC(am3$model), AIC(am0.2$model), AIC(am3.2$model),AIC(am3.3$model))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# save it
gtsave(selection_table_ab,"SimpleLUIModel_Selection_Abund.png", path = outDir)

# save model output tables for use in supplementary information 
# use function from sjPlot library to save neat versions of model output table
tab_model(am3.3$model, transform = NULL, file = paste0(outDir,"Output_table_abund.html"))
summary(am3.3$model) # check the table against the outputs
R2GLMER(am3.3$model) # check the R2 values 
# $conditional
# [1] 0.7623785
# 
# $marginal
# [1] 0.06004844

tab_model(sm3.3$model, transform = NULL, file = paste0(outDir,"Output_table_rich.html"))
summary(sm3.3$model) # check the table against the outputs
R2GLMER(sm3.3$model) # check the R2 values 
# $conditional
# [1] 0.7106989
# 
# $marginal
# [1] 0.05387385


# save model output tables for use in supplementary information 
# use function from sjPlot library to save neat versions of model output table
tab_model(am3$model, transform = NULL, file = paste0(outDir,"Output_table_abund_noOrder.html"))
summary(am3$model) # check the table against the outputs
R2GLMER(am3$model) # check the R2 values 
# $conditional
# [1] 0.6838617
# 
# $marginal
# [1] 0.007807586

tab_model(sm3$model, transform = NULL, file = paste0(outDir,"Output_table_rich_noOrder.html"))
summary(sm3$model) # check the table against the outputs
R2GLMER(sm3$model) # check the R2 values 
# $conditional
# [1] 0.6700185
# 
# $marginal
# [1] 0.007716608


##%######################################################%##
#                                                          #
####           Richness and abundance plots             ####
#                          Global                          #
#                                                          #
##%######################################################%##


####  5. Species Richness Plot ####
richness_metric <- predict_effects(iterations = 1000,
                                   model = sm3.3$model,
                                   model_data = model_data_sr,
                                   response_variable = "Species_richness",
                                   fixed_number = 2,
                                   fixed_column = c("Order", "LUI"),
                                   factor_number_1 = 5,
                                   factor_number_2 = 4,
                                   neg_binom = FALSE)

# rename prediction data frame and drop "Species_richness" column
result.sr <- fin_conf
result.sr <- dplyr::select(result.sr,-c(Species_richness))

richness_metric

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_sr$LUI), 5)[1:20]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric)

# plot species richness alone
richness <- richness_metric + xlab(NULL) + 
  labs(y ="Species Richness diff. (%)", x = "Order") +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_blank(), 
        legend.position = "right",
        legend.box = "vertical",
        legend.text = element_text(size = 7),
        legend.title = element_blank())

# save plot (pdf)
ggsave(filename = paste0(outDir, "SimpleLUI_Rich.pdf"), plot = last_plot(), width = 250, height = 100, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("SimpleLUI_Rich.jpeg", device ="jpeg", path = outDir, width=25, height=10, units="cm", dpi = 350)

####  6. Abundance Plot ####

abundance_metric <- predict_effects(iterations = 1000,
                                    model = am3.3$model,
                                    model_data = model_data_ab,
                                    response_variable = "LogAbund",
                                    fixed_number = 2,
                                    fixed_column = c("Order", "LUI"),
                                    factor_number_1 = 5,
                                    factor_number_2 = 4,
                                    neg_binom = FALSE)

# rename prediction data frame and drop "Abundance" column
result.ab <- fin_conf
result.ab <- dplyr::select(result.ab,-c(LogAbund))


abundance_metric

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab$LUI), 5)[1:20]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric)

# plot abundance alone
abundance <- abundance_metric + xlab(NULL) + 
  labs(y ="Total abundance diff. (%)", x = "Order") +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_blank(), 
        legend.position = "right",
        legend.box = "vertical",
        legend.text = element_text(size = 7),
        legend.title = element_blank())

# save plot (pdf)
ggsave(filename = paste0(outDir, "SimpleLUI_Abund.pdf"), plot = last_plot(), width = 250, height = 100, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("SimpleLUI_Abund.jpeg", device ="jpeg", path = outDir, width=25, height=10, units="cm", dpi = 350)

####  7. Plot Together  ####

richness <- richness_metric + 
  labs(y ="Species richness diff. (%)", x = NULL) +
  guides(scale = "none") +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "none")+
  ggtitle("a")


abundance <- abundance_metric +
  labs(y ="Total abundance diff. (%)", x = "Order") +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "none")+
  ggtitle("b")


# get the legend
legend <- get_legend(
  abundance +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right",
          legend.box = "vertical", 
          legend.text = element_text(size = 7), 
          legend.title = element_blank())
)

# plot together
plot_figure <- cowplot::plot_grid(richness, abundance, ncol = 1, rel_heights = c(1,1))
legend <- cowplot::plot_grid(NULL,legend,NULL, ncol = 1, rel_heights = c(0.5,1,0.5))
plot_figure <- cowplot::plot_grid(plot_figure, legend, ncol = 2, rel_widths = c(1,0.2))

# save plot (pdf)
ggsave(filename = paste0(outDir, "SimpleLUI.pdf"), plot = last_plot(), width = 250, height = 200, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("SimpleLUI.jpeg", device ="jpeg", path = outDir, width=20, height=15, units="cm", dpi = 350)

#### 8. Table of predicted values ####

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)
all_res$measure <- c(rep("ab", 20), rep("sr", 20))

# save as .csv
write.csv(all_res, file = paste0(predsDir,"LUI_predictions.csv"))

# save as table
LUI_predictions <- all_res %>% gt()
gtsave(LUI_predictions,"LUI_predictions.png", path = predsDir)

# # save as html, if wanted
# gtsave(LUI_predictions, paste0(predsDir, "LUI_predictions.html"))


############################################################
#                                                          #
#             Plots for models excluding Order             #
#                                                          #
############################################################

# Kyra's function above does not seem to work under these circumstances
## need to fix



# 
# #### Species Richness Plot ####
# richness_metric <- predict_effects(iterations = 1000,
#                                    model = sm3$model,
#                                    model_data = model_data_sr,
#                                    response_variable = "Species_richness",
#                                    fixed_number = 1,
#                                    fixed_column = "LUI",
#                                    factor_number_1 = 4,
#                                    neg_binom = FALSE)
# 
# # rename prediction data frame and drop "Species_richness" column
# result.sr <- fin_conf
# result.sr <- dplyr::select(result.sr,-c(Species_richness))
# 
# richness_metric
# 
# model_data <- function(model_plot){
#   ggplot_build(model_plot)$data[[2]] %>%
#     dplyr::select(y) %>%
#     cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
#     mutate(LUI = rep(levels(model_data_sr$LUI), 5)[1:20]) %>%
#     dplyr::select(LUI, y, ymin, ymax)
# }
# 
# # richness data
# model_data(richness_metric)
# 
# # plot species richness alone
# richness <- richness_metric + xlab(NULL) + 
#   labs(y ="Species Richness diff. (%)", x = "Order") +
#   scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
#   theme(axis.title = element_text(size = 8),
#         axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
#         axis.text.y = element_text(size = 7),
#         axis.ticks = element_blank(), 
#         legend.position = "right",
#         legend.box = "vertical",
#         legend.text = element_text(size = 7),
#         legend.title = element_blank())
# 
# # save plot (pdf)
# ggsave(filename = paste0(outDir, "SimpleLUI_Rich.pdf"), plot = last_plot(), width = 250, height = 100, units = "mm", dpi = 300)
# 
# # save plot (jpeg)
# ggsave("SimpleLUI_Rich.jpeg", device ="jpeg", path = outDir, width=25, height=10, units="cm", dpi = 350)
# 
# ####  6. Abundance Plot ####
# 
# abundance_metric <- predict_effects(iterations = 1000,
#                                     model = am3.3$model,
#                                     model_data = model_data_ab,
#                                     response_variable = "LogAbund",
#                                     fixed_number = 2,
#                                     fixed_column = c("Order", "LUI"),
#                                     factor_number_1 = 5,
#                                     factor_number_2 = 4,
#                                     neg_binom = FALSE)
# 
# # rename prediction data frame and drop "Abundance" column
# result.ab <- fin_conf
# result.ab <- dplyr::select(result.ab,-c(LogAbund))
# 
# 
# abundance_metric
# 
# model_data <- function(model_plot){
#   ggplot_build(model_plot)$data[[2]] %>%
#     dplyr::select(y) %>%
#     cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
#     mutate(LUI = rep(levels(model_data_ab$LUI), 5)[1:20]) %>%
#     dplyr::select(LUI, y, ymin, ymax)
# }
# 
# # abundance data
# model_data(abundance_metric)
# 
# # plot abundance alone
# abundance <- abundance_metric + xlab(NULL) + 
#   labs(y ="Total abundance diff. (%)", x = "Order") +
#   scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
#   theme(axis.title = element_text(size = 8),
#         axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
#         axis.text.y = element_text(size = 7),
#         axis.ticks = element_blank(), 
#         legend.position = "right",
#         legend.box = "vertical",
#         legend.text = element_text(size = 7),
#         legend.title = element_blank())
# 
# # save plot (pdf)
# ggsave(filename = paste0(outDir, "SimpleLUI_Abund.pdf"), plot = last_plot(), width = 250, height = 100, units = "mm", dpi = 300)
# 
# # save plot (jpeg)
# ggsave("SimpleLUI_Abund.jpeg", device ="jpeg", path = outDir, width=25, height=10, units="cm", dpi = 350)
# 
# ####  7. Plot Together  ####
# 
# richness <- richness_metric + 
#   labs(y ="Species richness diff. (%)", x = NULL) +
#   guides(scale = "none") +
#   scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
#         axis.text.y = element_text(size = 7),
#         legend.position = "none")+
#   ggtitle("a")
# 
# 
# abundance <- abundance_metric +
#   labs(y ="Total abundance diff. (%)", x = "Order") +
#   scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
#   theme(axis.title = element_text(size = 8),
#         axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
#         axis.text.y = element_text(size = 7),
#         legend.position = "none")+
#   ggtitle("b")
# 
# 
# # get the legend
# legend <- get_legend(
#   abundance +
#     guides(color = guide_legend(ncol = 1)) +
#     theme(legend.position = "right",
#           legend.box = "vertical", 
#           legend.text = element_text(size = 7), 
#           legend.title = element_blank())
# )
# 
# # plot together
# plot_figure <- cowplot::plot_grid(richness, abundance, ncol = 1, rel_heights = c(1,1))
# legend <- cowplot::plot_grid(NULL,legend,NULL, ncol = 1, rel_heights = c(0.5,1,0.5))
# plot_figure <- cowplot::plot_grid(plot_figure, legend, ncol = 2, rel_widths = c(1,0.2))
# 
# # save plot (pdf)
# ggsave(filename = paste0(outDir, "SimpleLUI.pdf"), plot = last_plot(), width = 250, height = 200, units = "mm", dpi = 300)
# 
# # save plot (jpeg)
# ggsave("SimpleLUI.jpeg", device ="jpeg", path = outDir, width=20, height=15, units="cm", dpi = 350)



# create dataframe for values to predict response to
nd <- data.frame(LUI=factor(c("Primary vegetation","Secondary vegetation",
                              "Agriculture_Low","Agriculture_High"), levels = levels(sm3$data$LUI)),
                 Species_richness=0,
                 LogAbund=0)

## species richness predictions ##
s.preds <- StatisticalModels::PredictGLMERRandIter(model = sm3$model, data = nd)

s.preds <- exp(s.preds)

# convert to percentage difference from primary vegetation
s.preds <- sweep(x = s.preds, MARGIN = 2, STATS = s.preds[1,], FUN = '/')

# get quantiles
s.preds.median <- ((apply(X = s.preds,MARGIN = 1,FUN = median))*100)-100
s.preds.upper <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
s.preds.lower <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



## abundance predictions ##
a.preds <- PredictGLMERRandIter(model = am3$model,data = nd)

a.preds <- exp(a.preds)-1

# convert to percentage difference from primary vegetation
a.preds <- sweep(x = a.preds,MARGIN = 2,STATS = a.preds[1,],FUN = '/')

# get quantiles
a.preds.median <- ((apply(X = a.preds,MARGIN = 1,FUN = median))*100)-100
a.preds.upper <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
a.preds.lower <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


# combine data into one table for plotting
abun_res <- as.data.frame(cbind(a.preds.median, a.preds.lower, a.preds.upper))
rich_res <- as.data.frame(cbind(s.preds.median, s.preds.lower, s.preds.upper))
colnames(abun_res) <- c("median", "lower", "upper")
colnames(rich_res) <- c("median", "lower", "upper")
abun_res$metric <- "abun"
rich_res$metric <- "rich"
abun_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))
rich_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))

abun_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA
rich_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA


# point points and error bars

p1 <- ggplot(data = abun_res) +
  geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
  geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  xlab("") +
  ylab("Change in total abundance (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
  ylim(c(-100, 0)) +
  theme(legend.position = "none", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2)) +
  ggtitle("b")


p2 <- ggplot(data = rich_res) +
  geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
  geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  xlab("") +
  ylab("Change in species richness (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
  ylim(c(-100, 0)) +
  theme(legend.position = "none", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2)) +
  ggtitle("a")


p3 <- cowplot::plot_grid(p2, p1)

ggsave(filename = paste0(outDir, "FIG_No_Order_LandUseOnly.pdf"), plot = last_plot(), width = 160, height = 90, units = "mm", dpi = 300)







