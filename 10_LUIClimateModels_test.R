
############################################################
#                                                          #
#    Additional tests: Compare AICs of models    #
#                                                          #
############################################################

# AIC values are used to compare models with and without Order as a fixed effect 

## Get set up ##

# directories 
inDir <- "5_RunLUIClimateModels/"
inDirTrop <- "6_TropicalModels/"
outDir <- "10_Additional_Tests/Model_selection/"
outDirTrop <- "10_Additional_Tests/Model_selection/Tropical/"
if(!dir.exists(inDir)) dir.create(inDir)
if(!dir.exists(inDirTrop)) dir.create(inDirTrop)
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(outDirTrop)) dir.create(outDirTrop)


sink(paste0(outDir,"log_LUI_ClimateModels_test.txt"))

t.start <- Sys.time()

print(t.start)

# load libraries
packages_model <- c("devtools","StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")

# read in the predicts climate data and models
predictsSites <- readRDS(paste0(inDir,"PREDICTSSitesClimate_Data.rds"))
load(paste0(inDir, "MeanAnomalyModelAbund.rdata"))
load(paste0(inDir, "MeanAnomalyModelRich.rdata"))
load(paste0(inDir, "MaxAnomalyModelAbund.rdata"))
load(paste0(inDir, "MaxAnomalyModelRich.rdata"))

## Model Selection - GLOBAL ##

# 1. Abundance, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


# MeanAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
#                                fixedStruct = "Order * LUI * StdTmeanAnomalyRS",
#                                randomStruct = "(1|SS)+(1|SSB)",
#                                saveVars = c("SSBS"))

MeanAnomalyModelAbund_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "LUI * StdTmeanAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))
# take a look at the AICs
AIC_MeanAbund<-print(AIC(MeanAnomalyModelAbund_test$model,MeanAnomalyModelAbund$model))

# save the model output
save(MeanAnomalyModelAbund_test, file = paste0(outDir, "MeanAnomalyModelAbund_test.rdata"))

# 2. Richness, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]

# MeanAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
#                               fixedStruct = "Order * LUI * StdTmeanAnomalyRS",
#                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
#                               saveVars = c("SSBS"))

MeanAnomalyModelRich_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "LUI * StdTmeanAnomalyRS",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS"))

# take a look at the AICs
AIC_MeanRich<-print(AIC(MeanAnomalyModelRich_test$model,MeanAnomalyModelRich$model))

# save the model output
save(MeanAnomalyModelRich_test, file = paste0(outDir, "MeanAnomalyModelRich_test.rdata"))

# 3. Abundance, max anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


# MaxAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
#                               fixedStruct = "Order * LUI * StdTmaxAnomalyRS",
#                               randomStruct = "(1|SS)+(1|SSB)",
#                               saveVars = c("SSBS"))

MaxAnomalyModelAbund_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                              fixedStruct = "LUI * StdTmaxAnomalyRS",
                              randomStruct = "(1|SS)+(1|SSB)",
                              saveVars = c("SSBS"))

# take a look at the AICs
AIC_MaxAbund<-print(AIC(MaxAnomalyModelAbund_test$model,MaxAnomalyModelAbund$model))

# save the model output
save(MaxAnomalyModelAbund_test, file = paste0(outDir, "MaxAnomalyModelAbund_test.rdata"))

# 4. Richness, max anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS),]

# MaxAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
#                              fixedStruct = "Order * LUI * StdTmaxAnomalyRS",
#                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
#                              saveVars = c("SSBS"))

MaxAnomalyModelRich_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                             fixedStruct = "LUI * StdTmaxAnomalyRS",
                             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                             saveVars = c("SSBS"))

# take a look at the AICs
AIC_MaxRich<-print(AIC(MaxAnomalyModelRich_test$model, MaxAnomalyModelRich$model))

# save the model output
save(MaxAnomalyModelRich_test, file = paste0(outDir, "MaxAnomalyModelRich_test.rdata"))

# put together and save .csv file
AICs_test <- rbind(AIC_MeanAbund,AIC_MaxAbund,AIC_MeanRich,AIC_MaxRich)

write.csv(AICs_test, file = paste0(outDir,"AICs_test.csv"))

#### save model output tables for use in supplementary information ####
# use function from sjPlot library to save neat versions of model output table
# conditional: the conditional R2 value, i.e. the variance explained by fixed and random effects 
# marginal: the marginal R2 value, i.e. the variance explained by the fixed effects

tab_model(MeanAnomalyModelAbund_test$model, transform = NULL, file = paste0(outDir,"Tables/AbunMeanAnom_test_output_table.html"))
summary(MeanAnomalyModelAbund_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_test$model) # check the R2 values 
# $conditional
# [1] 0.3417518
# 
# $marginal
# [1] 0.03110523

tab_model(MeanAnomalyModelRich_test$model, transform = NULL, file = paste0(outDir,"Tables/RichMeanAnom_test_output_table.html"))
summary(MeanAnomalyModelRich_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_test$model) # check the R2 values
# $conditional
# [1] 0.6715986
# 
# $marginal
# [1] 0.009126403

tab_model(MaxAnomalyModelAbund_test$model, transform = NULL, file = paste0(outDir,"Tables/AbunMaxAnom_test_output_table.html"))
summary(MaxAnomalyModelAbund_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_test$model) # check the R2 values 
# $conditional
# [1] 0.337773
# 
# $marginal
# [1] 0.025435

tab_model(MaxAnomalyModelRich_test$model, transform = NULL, file = paste0(outDir,"Tables/RichMaxAnom_test_output_table.html"))
summary(MaxAnomalyModelRich_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_test$model) # check the R2 values 
# $conditional
# [1] 0.6688694
# 
# $marginal
# [1] 0.008126281

# table of AICs
# species richness and abundance together

# load AIC.csv file
AICs <- read.csv("C:/Users/Kyra/Documents/GLITRS/Code/10_Additional_Tests/AICs_test.csv", header=TRUE, stringsAsFactors=FALSE)
# AIC values pulled from CSV file

# selection_table <- data.frame("Response" = c(rep("Species richness", 4),
#                                              rep("Total abundance", 4)),
#                               "Climate Anomaly" = c("Mean","Mean","Maximum","Maximum",
#                                                     "Mean","Mean","Maximum","Maximum"),
#                               "Model" = c("Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Total_abundance ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)"),
#                               "AIC" = c(AIC(MeanAnomalyModelRich_test$model), AIC(MeanAnomalyModelRich$model), AIC(MaxAnomalyModelRich_test$model), AIC(MaxAnomalyModelRich$model),  
#                                         AIC(MeanAnomalyModelAbund_test$model), AIC(MeanAnomalyModelAbund$model), AIC(MaxAnomalyModelAbund_test$model), AIC(MaxAnomalyModelAbund$model))) %>%
#   group_by(Response) %>%                              
#   mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
#   ungroup() %>%
#   gt()

# add model descriptions
AICs$Model <- c("Abundance ~ LUI * STA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * STA + (1|SS) + (1|SSB)",
                "Abundance ~ LUI * SMTA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * SMTA + (1|SS) + (1|SSB)",
                "Species richness ~ LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ LUI * SMTA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * SMTA + (1|SS) + (1|SSB) + (1|SSBS)")


# add climate anomaly
AICs$Anomaly <- c(rep("Mean",2),rep("Maximum",2),rep("Mean",2),rep("Maximum",2))

# add climate anomaly
AICs$Response <- c(rep("Abundance",4),rep("Species_richness",4))

AICs <- as.data.frame(AICs)

# drop columns 'X' and 'df' and re-order columns
AICs <- subset (AICs,select = c(Response,Anomaly,Model,AIC))

# make gt table of model selection
# AIC_select <- AICs %>% 
#   group_by(Response,Anomaly) %>%                              
#   mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
#   ungroup() %>%
#   select(Anomaly,Model,AIC,deltaAIC) %>%
#   gt(rowname_col = "Model") %>%
#   tab_row_group(
#     label = "Species Richness",
#     rows = starts_with("Species richness")
#   ) %>%
#   tab_row_group(
#     label = "Abundance",
#     rows = starts_with("Abundance")
#   ) %>%
#   tab_stubhead(label = "Models")

AIC_select <- AICs %>% 
  group_by(Response,Anomaly) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  select(Model,Anomaly,AIC,deltaAIC) %>%
  gt(rowname_col = "Model") %>%
  tab_row_group(
    label = "Species Richness",
    rows = starts_with("Species richness")
  ) %>%
  tab_row_group(
    label = "Abundance",
    rows = starts_with("Abundance")
  ) %>% 
  cols_align(
    align = "center",
    columns = c(Model, Anomaly, AIC, deltaAIC)
  )%>%
  tab_stubhead(label = "Models") %>%
  cols_label(
    Model = md ("Models"),
    Anomaly = md("Temperature Anomaly"),
    AIC = md("AIC"),
    deltaAIC = md("deltaAIC")
  ) 

# save
# sometimes gtsave requires the complete file pathway
gtsave(AIC_select,"C:/Users/Kyra/Documents/GitHub/LUCC_ByOrder_paper/10_Additional_Tests/Model_selection/AIC_ModelSelection.png")

## Model Selection - TROPICAL/NONTROPICAL ##

# read in data and model
nontrop <- readRDS(paste0(inDirTrop,"nontrop.rds"))
trop <- readRDS(paste0(inDirTrop,"trop.rds"))
load(paste0(inDirTrop, "MeanAnomalyModelAbund_nontrop.rdata"))
load(paste0(inDirTrop, "MeanAnomalyModelAbund_trop.rdata"))
load(paste0(inDirTrop, "MeanAnomalyModelRich_nontrop.rdata"))
load(paste0(inDirTrop, "MeanAnomalyModelRich_trop.rdata"))
load(paste0(inDirTrop, "MaxAnomalyModelAbund_nontrop.rdata"))
load(paste0(inDirTrop, "MaxAnomalyModelAbund_trop.rdata"))
load(paste0(inDirTrop, "MaxAnomalyModelRich_nontrop.rdata"))
load(paste0(inDirTrop, "MaxAnomalyModelRich_trop.rdata"))

# 5a. Abundance, Mean Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund_nontrop_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                       fixedStruct = "LUI * StdTmeanAnomalyRS",
                                       randomStruct = "(1|SS)+(1|SSB)",
                                       saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund_nontrop_test$model)

# take a look at the AICs
AIC_MeanAbund_nontrop<-print(AIC(MeanAnomalyModelAbund_nontrop_test$model,MeanAnomalyModelAbund_nontrop$model))

# save the model output
save(MeanAnomalyModelAbund_nontrop_test, file = paste0(outDirTrop, "MeanAnomalyModelAbund_nontrop_test.rdata"))

# 6a. Richness, Mean Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich_nontrop_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                      fixedStruct = "LUI * StdTmeanAnomalyRS",
                                      randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                      saveVars = c("SSBS"))
# get summary
summary(MeanAnomalyModelRich_nontrop_test$model)

# take a look at the AICs
AIC_MeanRich_nontrop<-print(AIC(MeanAnomalyModelRich_nontrop_test$model,MeanAnomalyModelRich_nontrop$model))

# save the model output
save(MeanAnomalyModelRich_nontrop_test, file = paste0(outDirTrop, "MeanAnomalyModelRich_nontrop_test.rdata"))

# 7a. Abundance, Max Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund_nontrop_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                      fixedStruct = "LUI * StdTmaxAnomalyRS",
                                      randomStruct = "(1|SS)+(1|SSB)",
                                      saveVars = c("SSBS"))
# get summary
summary(MaxAnomalyModelAbund_nontrop_test$model)

# take a look at the AICs
AIC_MaxAbund_nontrop<-print(AIC(MaxAnomalyModelAbund_nontrop_test$model,MaxAnomalyModelAbund_nontrop$model))

# save the model output
save(MaxAnomalyModelAbund_nontrop_test, file = paste0(outDirTrop, "MaxAnomalyModelAbund_nontrop_test.rdata"))

# 8a. Richness, Max Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$StdTmeanAnomalyRS),]

MaxAnomalyModelRich_nontrop_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                     fixedStruct = "LUI * StdTmaxAnomalyRS",
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    saveVars = c("SSBS"))
# get summary
summary(MaxAnomalyModelRich_nontrop_test$model)

# take a look at the AICs
AIC_MaxRich_nontrop<-print(AIC(MaxAnomalyModelRich_nontrop_test$model,MaxAnomalyModelRich_nontrop$model))

# save the model output
save(MaxAnomalyModelRich_nontrop_test, file = paste0(outDirTrop, "MaxAnomalyModelRich_nontrop_test.rdata"))

# 5b. Abundance, mean anomaly, Tropical

model_data <- trop[!is.na(trop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund_trop_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                    fixedStruct = "LUI * StdTmeanAnomalyRS",
                                    randomStruct = "(1|SS)+(1|SSB)",
                                    saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund_trop_test$model)

# take a look at the AICs
AIC_MeanAbund_trop<-print(AIC(MeanAnomalyModelAbund_trop_test$model,MeanAnomalyModelAbund_trop$model))

# save the model output
save(MeanAnomalyModelAbund_trop_test, file = paste0(outDirTrop, "MeanAnomalyModelAbund_trop_test.rdata"))

# 6b. Richness, mean anomaly, Tropical

model_data <- trop[!is.na(trop$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich_trop_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                   fixedStruct = "LUI * StdTmeanAnomalyRS",
                                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                   saveVars = c("SSBS"))
# get summary
summary(MeanAnomalyModelRich_trop_test$model)

# take a look at the AICs
AIC_MeanRich_trop<-print(AIC(MeanAnomalyModelRich_trop_test$model,MeanAnomalyModelRich_trop$model))

# save the model output
save(MeanAnomalyModelRich_trop_test, file = paste0(outDirTrop, "MeanAnomalyModelRich_trop_test.rdata"))

# 7b. Abundance, Max Anomaly, Tropical

model_data <- trop[!is.na(trop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund_trop_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                   fixedStruct = "LUI * StdTmaxAnomalyRS",
                                   randomStruct = "(1|SS)+(1|SSB)",
                                   saveVars = c("SSBS"))
# get summary
summary(MaxAnomalyModelAbund_trop_test$model)

# take a look at the AICs
AIC_MaxAbund_trop<-print(AIC(MaxAnomalyModelAbund_trop_test$model,MaxAnomalyModelAbund_trop$model))

# save the model output
save(MaxAnomalyModelAbund_trop_test, file = paste0(outDirTrop, "MaxAnomalyModelAbund_trop_test.rdata"))

# 8b. Richness, Max Anomaly, Tropical

model_data <- trop[!is.na(trop$StdTmeanAnomalyRS),]

MaxAnomalyModelRich_trop_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                  fixedStruct = "LUI * StdTmaxAnomalyRS",
                                  randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                  saveVars = c("SSBS"))
# get summary
summary(MaxAnomalyModelRich_trop_test$model)

# take a look at the AICs
AIC_MaxRich_trop<-print(AIC(MaxAnomalyModelRich_trop_test$model,MaxAnomalyModelRich_trop$model))

# save the model output
save(MaxAnomalyModelRich_trop_test, file = paste0(outDirTrop, "MaxAnomalyModelRich_trop_test.rdata"))


# put together AICs and save .csv file
AICs_trop_test <- rbind(AIC_MeanAbund_nontrop, AIC_MeanAbund_trop,
                        AIC_MeanRich_nontrop, AIC_MeanRich_trop,
                        AIC_MaxAbund_nontrop, AIC_MaxAbund_trop,
                        AIC_MaxRich_nontrop, AIC_MaxRich_trop)

# save
write.csv(AICs_trop_test, file = paste0(outDirTrop,"AICs_trop_test.csv"))

#### save model output tables for use in supplementary information ####
# use function from sjPlot library to save neat versions of model output table
# conditional: the conditional R2 value, i.e. the variance explained by fixed and random effects 
# marginal: the marginal R2 value, i.e. the variance explained by the fixed effects

tab_model(MeanAnomalyModelAbund_nontrop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/AbunMeanAnom_nontrop_test_output_table.html"))
summary(MeanAnomalyModelAbund_nontrop_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_nontrop_test$model) # check the R2 values 
# $conditional
# [1] 0.4041423
# 
# $marginal
# [1] 0.00751681

tab_model(MeanAnomalyModelRich_nontrop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/RichMeanAnom_nontrop_test_output_table.html"))
summary(MeanAnomalyModelRich_nontrop_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_nontrop_test$model) # check the R2 values
# $conditional
# [1] 0.687376
# 
# $marginal
# [1] 0.01107587

tab_model(MaxAnomalyModelAbund_nontrop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/AbunMaxAnom_nontrop_test_output_table.html"))
summary(MaxAnomalyModelAbund_nontrop_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_nontrop_test$model) # check the R2 values 
# $conditional
# [1] 0.4011494
# 
# $marginal
# [1] 0.008403599

tab_model(MaxAnomalyModelRich_nontrop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/RichMaxAnom_nontrop_test_output_table.html"))
summary(MaxAnomalyModelRich_nontrop_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_nontrop_test$model) # check the R2 values 
# $conditional
# [1] 0.6759265
# 
# $marginal
# [1] 0.005548422

tab_model(MeanAnomalyModelAbund_trop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/AbunMeanAnom_trop_test_output_table.html"))
summary(MeanAnomalyModelAbund_trop_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_trop_test$model) # check the R2 values 
# $conditional
# [1] 0.3225042
# 
# $marginal
# [1] 0.07819674

tab_model(MeanAnomalyModelRich_trop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/RichMeanAnom_trop_test_output_table.html"))
summary(MeanAnomalyModelRich_trop_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_trop_test$model) # check the R2 values
# $conditional
# [1] 0.628909
# 
# $marginal
# [1] 0.02340666

tab_model(MaxAnomalyModelAbund_trop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/AbunMaxAnom_trop_test_output_table.html"))
summary(MaxAnomalyModelAbund_trop_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_trop_test$model) # check the R2 values 
# $conditional
# [1] 0.3219352
# 
# $marginal
# [1] 0.06952179

tab_model(MaxAnomalyModelRich_trop_test$model, transform = NULL, file = paste0(outDirTrop,"Tables/RichMaxAnom_trop_test_output_table.html"))
summary(MaxAnomalyModelRich_trop_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_trop_test$model) # check the R2 values 
# $conditional
# [1] 0.6291205
# 
# $marginal
# [1] 0.01835229

# table of AICs
# species richness and abundance together

# load AIC.csv file, if necessary
# AICs_trop_test <- read.csv(outDirTrop, "AICs_trop_test.csv", header=TRUE, stringsAsFactors=FALSE)

# selection_table <- data.frame("Response" = c(rep("Species richness", 4),
#                                              rep("Total abundance", 4)),
#                               "Climate Anomaly" = c("Mean","Mean","Maximum","Maximum",
#                                                     "Mean","Mean","Maximum","Maximum"),
#                               "Model" = c("Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Total_abundance ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)"),
#                               "AIC" = c(AIC(MeanAnomalyModelRich_test$model), AIC(MeanAnomalyModelRich$model), AIC(MaxAnomalyModelRich_test$model), AIC(MaxAnomalyModelRich$model),  
#                                         AIC(MeanAnomalyModelAbund_test$model), AIC(MeanAnomalyModelAbund$model), AIC(MaxAnomalyModelAbund_test$model), AIC(MaxAnomalyModelAbund$model))) %>%
#   group_by(Response) %>%                              
#   mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
#   ungroup() %>%
#   gt()

# add model descriptions
AICs_trop_test$Model <- c("Abundance ~ LUI * STA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * STA + (1|SS) + (1|SSB)",
                "Abundance ~ LUI * STA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * STA + (1|SS) + (1|SSB)",
                "Species richness ~ LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Abundance ~ LUI * SMTA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * SMTA + (1|SS) + (1|SSB)",
                "Abundance ~ LUI * SMTA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * SMTA + (1|SS) + (1|SSB)",
                "Species richness ~ LUI * SMTA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * SMTA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)")


# add climate anomaly
AICs_trop_test$Anomaly <- c(rep("Mean",8),rep("Maximum",8))

# add diversity variable
AICs_trop_test$Response <- c(rep("Abundance",4),
                             rep("Species_richness",4),
                             rep("Abundance",4),
                             rep("Species_richness",4))

# add realm variable
AICs_trop_test$Realm <- c(rep("NonTropical",2),rep("Tropical",2),
                          rep("NonTropical",2),rep("Tropical",2),
                          rep("NonTropical",2),rep("Tropical",2),
                          rep("NonTropical",2),rep("Tropical",2))

AICs_trop_test <- as.data.frame(AICs_trop_test)

# drop columns 'X' and 'df' and re-order columns
AICs_trop_test <- subset (AICs_trop_test,select = c(Response,Anomaly,Model,AIC))

AIC_select <- AICs_trop_test %>% 
  group_by(Response,Anomaly,Realm) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  select(Model,Realm,Anomaly,AIC,deltaAIC) %>%
  gt(rowname_col = "Model") %>%
  tab_row_group(
    label = "Species Richness",
    rows = starts_with("Species richness")
  ) %>%
  tab_row_group(
    label = "Abundance",
    rows = starts_with("Abundance")
  ) %>% 
  cols_align(
    align = "center",
    columns = c(Model, Realm, Anomaly, AIC, deltaAIC)
  )%>%
  tab_stubhead(label = "Models") %>%
  cols_label(
    Model = md ("Models"),
    Anomaly = md("Temperature Anomaly"),
    AIC = md("AIC"),
    deltaAIC = md("deltaAIC")
  ) 

# save
# sometimes gtsave requires the complete file pathway
gtsave(AIC_select,"C:/Users/Kyra/Documents/GitHub/LUCC_ByOrder_paper/10_Additional_Tests/Model_selection/Tropical/AIC_ModelSelection.png")

# make a model selection table for only MeanAnomaly models

AIC_select_STA <- AICs_trop_test %>% 
  filter(Anomaly == "Mean") %>%
  group_by(Response,Realm) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  select(Model,Realm,AIC,deltaAIC) %>%
  gt(rowname_col = "Model") %>%
  tab_row_group(
    label = "Species Richness",
    rows = starts_with("Species richness")
  ) %>%
  tab_row_group(
    label = "Abundance",
    rows = starts_with("Abundance")
  ) %>% 
  cols_align(
    align = "center",
    columns = c(Model, Realm, AIC, deltaAIC)
  )%>%
  tab_stubhead(label = "Models") %>%
  cols_label(
    Model = md ("Models"),
    AIC = md("AIC"),
    deltaAIC = md("deltaAIC")
  ) 
# save
# sometimes gtsave requires the complete file pathway
gtsave(AIC_select_STA,"C:/Users/Kyra/Documents/GitHub/LUCC_ByOrder_paper/10_Additional_Tests/Model_selection/Tropical/AIC_ModelSelection_STA.png")


t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
