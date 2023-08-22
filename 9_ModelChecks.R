##%######################################################%##
#                                                          #
####                    Model Checks                    ####
#                                                          #
##%######################################################%##

# This script carries out checks on the models including
# Q-Q plots, fitted vs residuals and and a test of spatial
# autocorrelations using Moran's I. 


rm(list = ls())

# directories
LUIDir <- "2_RunSimpleLUIModel/"
LUICCDir <- "5_RunLUIClimateModels/"
predictsDataDir <- "5_RunLUIClimateModels/"
outdir <- "9_ModelChecks/"
if(!dir.exists(outdir)) dir.create(outdir)

# load libraries
library(cowplot)
library(ggplot2)
library(StatisticalModels)
library(gridGraphics)

# read in model files
# land use only
LUIAbund <- readRDS(file = paste0(LUIDir,"am3.3.rds"))
LUIRich <- readRDS(file = paste0(LUIDir,"sm3.3.rds"))
load(file = paste0(LUIDir,"Abundance_landuse_model_noOrder.rdata")) # am3
load(file = paste0(LUIDir,"Richness_landuse_model_noOrder.rdata")) # sm3

# land use, STA and Order
load(paste0(LUICCDir, "MeanAnomalyModelAbund.rdata")) # MeanAnomalyModelAbund
load(paste0(LUICCDir, "MeanAnomalyModelRich.rdata")) # MeanAnomalyModelRich

# land use and STA only
load(paste0(LUICCDir, "MeanAnomalyModelAbund_noOrder.rdata")) # MeanAnomalyModelAbund2
load(paste0(LUICCDir, "MeanAnomalyModelRich_noOrder.rdata")) # MeanAnomalyModelRich2


# try to write a function to run these checks on each model output and save a pdf of all figs.

#mod_list <- c("LUIAbund", "LUIRich","MeanAnomalyModelAbund", "MeanAnomalyModelRich", "MeanAnomalyModelAbund2", "MeanAnomalyModelRich2")
mod_list1 <- c("LUIAbund", "LUIRich")
mod_list2 <- c("MeanAnomalyModelAbund", "MeanAnomalyModelRich", "MeanAnomalyModelAbund2", "MeanAnomalyModelRich2")

# load datasets
model_data_sr <- readRDS(file = paste0(LUIDir,"model_data_sr.rds"))
model_data_ab <- readRDS(file = paste0(LUIDir,"model_data_ab.rds"))

predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesClimate_Data.rds"))

#x <- mod_list1[1]

# for the Order*LUI models
for(x in mod_list1){
  
  # only do this one if not SR model
  # grepl() searches for matches to "Abun" in "x", the list of models
  if(grepl("Abun", x) == TRUE) {
    
    
    ## 1. Checking the fitted vs residuals relationship
    p1 <- plot(get(x)$model)
  }
  
  ## 2. Normality of Residuals
  pdf(NULL)
  dev.control(displaylist="enable")
  qqnorm(resid(get(x)$model), main = "")
  qqline(resid(get(x)$model))
  p2 <- recordPlot()
  invisible(dev.off())
  
  
  if(grepl("Abun", x) == TRUE) {
    
    
    predData <- model_data_ab[!is.na(model_data_ab$LogAbund), ]
    
    
    # 3. plot of observed vs fitted values
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predData$LogAbund,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    
    cowplot::plot_grid(p1,p2,p3,
                       labels = c("(a)", "(b)", "(c)"))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3)
    
    
  }else{
    
    
    # 4. plot of observed vs fitted values
    
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(model_data_sr$Species_richness,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    cowplot::plot_grid(p2,p3,
                       labels = c("A.", "B."))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf")) 
    
    rm(p2, p3)
    
  }
  
}

# for the Order*LUI*Climate Models

# loop through models and generate plots for checking assumptions etc
for(x in mod_list2){
  
  # only do this one if not SR model
  # grepl() searches for matches to "Abun" in "x", the list of models
  if(grepl("Abun", x) ==1) {
    
    
    ## 1. Checking the fitted vs residuals relationship
    p1 <- plot(get(x)$model)
  }
  
  ## 2. Normality of Residuals
  pdf(NULL)
  dev.control(displaylist="enable")
  qqnorm(resid(get(x)$model), main = "")
  qqline(resid(get(x)$model))
  p2 <- recordPlot()
  invisible(dev.off())
  
  
  if(grepl("Abun", x) == 1) {
    
    
    predData <- predictsSites[!is.na(predictsSites$LogAbund), ]
    
    
    # 3. plot of observed vs fitted values
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predData$LogAbund,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    
    cowplot::plot_grid(p1,p2,p3,
                       labels = c("A.", "B.", "C."))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3)
    
    
  }else{
    
    
    # 4. plot of observed vs fitted values
    
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predictsSites$Species_richness,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    cowplot::plot_grid(p2,p3,
                       labels = c("A.", "B."))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf")) 
    
    rm(p2, p3)
    
  }
  
}


############################################################
#                                                          #
#            Model checks alternative packages             #
#                                                          #
############################################################

## using the performance package
# https://easystats.github.io/performance/articles/check_model.html
library(performance)

check_model(LUIRich$model)
check_model(LUIAbund$model)

check_model(MeanAnomalyModelAbund$model)
check_model(MeanAnomalyModelAbund2$model)

check_model(MeanAnomalyModelRich$model)
check_model(MeanAnomalyModelRich2$model)

## using the DHARMa package
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#interpreting-residuals-and-recognizing-misspecification-problems
library(DHARMa)


#### LUI only models ####

# calculate residuals first to avoid rerunning within each function
simulationOutput_richLUI <- simulateResiduals(fittedModel = LUIRich$model, plot = F)
simulationOutput_abunLUI <- simulateResiduals(fittedModel = LUIAbund$model, plot = F)

simulationOutput_richLUI_NO <- simulateResiduals(fittedModel = sm3$model, plot = F)
simulationOutput_abunLUI_NO <- simulateResiduals(fittedModel = am3$model, plot = F)

# Q-Q plot and residuals vs predicted
plot(simulationOutput_richLUI)
plot(simulationOutput_abunLUI)

plot(simulationOutput_richLUI_NO)
plot(simulationOutput_abunLUI_NO)

# residuals vs predictor
plotResiduals(simulationOutput_richLUI, form = model_data_sr$LUI)
plotResiduals(simulationOutput_richLUI, form = model_data_sr$Order)
plotResiduals(simulationOutput_abunLUI, form = model_data_ab$LUI)
plotResiduals(simulationOutput_abunLUI, form = model_data_ab$Order)


plotResiduals(simulationOutput_richLUI_NO, form = model_data_sr$LUI)
plotResiduals(simulationOutput_abunLUI_NO, form = model_data_ab$LUI)

# test uniformity
testUniformity(simulationOutput_richLUI)
testUniformity(simulationOutput_abunLUI)

testUniformity(simulationOutput_richLUI_NO)
testUniformity(simulationOutput_abunLUI_NO)

# test outliers
testOutliers(simulationOutput_richLUI)
testOutliers(simulationOutput_abunLUI) # abundance has outliers but not clear where

testOutliers(simulationOutput_richLUI_NO)
testOutliers(simulationOutput_abunLUI_NO) # abundance has outliers but not clear where

# test dispersion
testDispersion(simulationOutput_richLUI) 
testDispersion(simulationOutput_abunLUI) 
testDispersion(simulationOutput_richLUI, type = "PearsonChisq", alternative = "greater") 
testDispersion(simulationOutput_abunLUI, type = "PearsonChisq", alternative = "greater") 

testDispersion(simulationOutput_richLUI_NO) 
testDispersion(simulationOutput_abunLUI_NO) 
testDispersion(simulationOutput_richLUI_NO, type = "PearsonChisq", alternative = "greater") 
testDispersion(simulationOutput_abunLUI_NO, type = "PearsonChisq", alternative = "greater") 

# test zero inflation
testZeroInflation(simulationOutput_richLUI)
testZeroInflation(simulationOutput_abunLUI)

testZeroInflation(simulationOutput_richLUI_NO)
testZeroInflation(simulationOutput_abunLUI_NO)

#### LUI and STA models ####

# calculate residuals first to avoid rerunning within each function
simulationOutput_richLUISTA <- simulateResiduals(fittedModel = MeanAnomalyModelRich$model, plot = F)
simulationOutput_abunLUISTA <- simulateResiduals(fittedModel = MeanAnomalyModelAbund$model, plot = F)


simulationOutput_richLUISTA_NO <- simulateResiduals(fittedModel = MeanAnomalyModelRich2$model, plot = F)
simulationOutput_abunLUISTA_NO <- simulateResiduals(fittedModel = MeanAnomalyModelAbund2$model, plot = F)

# Q-Q plot and residuals vs predicted
plot(simulationOutput_richLUISTA)
plot(simulationOutput_abunLUISTA)

plot(simulationOutput_richLUISTA_NO)
plot(simulationOutput_abunLUISTA_NO)

# residuals vs predictor
plotResiduals(simulationOutput_richLUISTA, form = MeanAnomalyModelRich$data$LUI)
plotResiduals(simulationOutput_richLUISTA, form = MeanAnomalyModelRich$data$Order)
plotResiduals(simulationOutput_richLUISTA, form = MeanAnomalyModelRich$data$StdTmeanAnomalyRS)

plotResiduals(simulationOutput_abunLUISTA, form = MeanAnomalyModelAbund$data$LUI)
plotResiduals(simulationOutput_abunLUISTA, form = MeanAnomalyModelAbund$data$Order)
plotResiduals(simulationOutput_abunLUISTA, form = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS)

plotResiduals(simulationOutput_richLUISTA_NO, form = MeanAnomalyModelRich2$data$LUI)
plotResiduals(simulationOutput_richLUISTA_NO, form = MeanAnomalyModelRich2$data$StdTmeanAnomalyRS)

plotResiduals(simulationOutput_abunLUISTA_NO, form = MeanAnomalyModelAbund2$data$LUI)
plotResiduals(simulationOutput_abunLUISTA_NO, form = MeanAnomalyModelAbund2$data$StdTmeanAnomalyRS)


# test uniformity
testUniformity(simulationOutput_richLUISTA)
testUniformity(simulationOutput_abunLUISTA)

testUniformity(simulationOutput_richLUISTA_NO)
testUniformity(simulationOutput_abunLUISTA_NO)

# test outliers
testOutliers(simulationOutput_richLUISTA)
testOutliers(simulationOutput_abunLUISTA) 

testOutliers(simulationOutput_richLUISTA_NO)
testOutliers(simulationOutput_abunLUISTA_NO)

# test dispersion
testDispersion(simulationOutput_richLUISTA) 
testDispersion(simulationOutput_abunLUISTA) 
testDispersion(simulationOutput_richLUISTA, type = "PearsonChisq", alternative = "greater") 
testDispersion(simulationOutput_abunLUISTA, type = "PearsonChisq", alternative = "greater") 

testDispersion(simulationOutput_richLUISTA_NO) 
testDispersion(simulationOutput_abunLUISTA_NO) 
testDispersion(simulationOutput_richLUISTA_NO, type = "PearsonChisq", alternative = "greater") 
testDispersion(simulationOutput_abunLUISTA_NO, type = "PearsonChisq", alternative = "greater") 

# test zero inflation
testZeroInflation(simulationOutput_richLUISTA)
testZeroInflation(simulationOutput_abunLUISTA)

testZeroInflation(simulationOutput_richLUISTA_NO)
testZeroInflation(simulationOutput_abunLUISTA_NO)



