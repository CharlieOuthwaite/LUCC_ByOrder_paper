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
# land use and Order
load(file = paste0(LUIDir,"Abundance_landuse_model.rdata")) # am3.3
load(file = paste0(LUIDir,"Richness_landuse_model.rdata")) # sm3.3

# land use only
load(file = paste0(LUIDir,"Abundance_landuse_model_noOrder.rdata")) # am3
load(file = paste0(LUIDir,"Richness_landuse_model_noOrder.rdata")) # sm3

# land use, STA and Order
load(paste0(LUICCDir, "MeanAnomalyModelAbund_outliersrm.rdata")) # MeanAnomalyModelAbund
load(paste0(LUICCDir, "MeanAnomalyModelRich_outliersrm.rdata")) # MeanAnomalyModelRich

# land use and STA only
load(paste0(LUICCDir, "MeanAnomalyModelAbund_noOrder_outliersrm.rdata")) # MeanAnomalyModelAbund2
load(paste0(LUICCDir, "MeanAnomalyModelRich_noOrder_outliersrm.rdata")) # MeanAnomalyModelRich2

# list models for loops
mod_list1 <- c("am3.3", "sm3.3", "am3", "sm3") # land use models
mod_list2 <- c("MeanAnomalyModelAbund", "MeanAnomalyModelRich", "MeanAnomalyModelAbund2", "MeanAnomalyModelRich2") #land use/sta models

# load datasets
model_data_sr <- readRDS(file = paste0(LUIDir,"model_data_sr.rds"))
model_data_ab <- readRDS(file = paste0(LUIDir,"model_data_ab.rds"))

predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesClimate_Data.rds"))

#x <- mod_list1[1]

# for the Order*LUI models
for(x in mod_list1){
  
  # only do this one if not SR model
  # grepl() searches for matches to "Abun" in "x", the list of models
  if(x %in% c("am3.3", "am3")) {
    
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
  
  
  if(x %in% c("am3.3", "am3")) {
    
    #predData <- model_data_ab[!is.na(model_data_ab$LogAbund), ]
    
    
    # 3. plot of observed vs fitted values
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(model_data_ab$LogAbund,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    
    
    ## 4. Check for spatial autocorrelation
    
    sa_test<-SpatialAutocorrelationTest(model=get(x), ranefGrouping = "SS")
    
    #summary(sa_test)
    
    # percentage of studies that show spatial autocorrelation?
    perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100
    
    
    sa_test_vals <- as.data.frame(sa_test$P)
    sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)
    
    label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")
    
    p4 <- ggplot(data = sa_test_vals ) +
      geom_histogram(aes(x =`sa_test$P`)) +
      geom_vline(xintercept = 0.05, col = "red") +
      geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
      theme_bw() +
      ylim(c(0, 100)) +
      xlab("P-value") +
      ylab("Frequency") +
      theme(panel.grid = element_blank(), 
            aspect.ratio = 1)
    
    
    cowplot::plot_grid(p1, p2, p3, p4,
                       labels = c("(a)", "(b)", "(c)", "(d)"))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3, p4)
    
    
  }else{
    
    
    # 4. plot of observed vs fitted values
    
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(model_data_sr$Species_richness,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    ## 4. Check for spatial autocorrelation
    
    sa_test<-SpatialAutocorrelationTest(model=get(x), ranefGrouping = "SS")
    
    #summary(sa_test)
    
    # percentage of studies that show spatial autocorrelation?
    perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100
    
    
    sa_test_vals <- as.data.frame(sa_test$P)
    sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)
    
    label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")
    
    p4 <- ggplot(data = sa_test_vals ) +
      geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
      geom_vline(xintercept = 0.05, col = "red") +
      geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
      theme_bw() +
      ylim(c(0, 100)) +
      xlab("P-value") +
      ylab("Frequency") +
      theme(panel.grid = element_blank(), 
            aspect.ratio = 1)
    
    cowplot::plot_grid(p2, p3, p4,
                       labels = c("(a)", "(b)", "(c)"))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p2, p3, p4)
    
  }
  
}




# for the Order*LUI*Climate Models

#x <- mod_list2[1]

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
    
    
    ## 4. Check for spatial autocorrelation
    
    sa_test<-SpatialAutocorrelationTest(model=get(x), ranefGrouping = "SS")
    
    #summary(sa_test)
    
    # percentage of studies that show spatial autocorrelation?
    perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100
    
    
    sa_test_vals <- as.data.frame(sa_test$P)
    sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)
    
    label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")
    
    p4 <- ggplot(data = sa_test_vals ) +
      geom_histogram(aes(x = `sa_test$P`)) +
      geom_vline(xintercept = 0.05, col = "red") +
      geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
      theme_bw() +
      ylim(c(0, 100)) +
      xlab("P-value") +
      ylab("Frequency") +
      theme(panel.grid = element_blank(), 
            aspect.ratio = 1)
    
    cowplot::plot_grid(p1, p2, p3, p4, 
                       labels = c("(a)", "(b)", "(c)", "(d)"))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3, p4)
    
    
    
  }else{
    
    
    # 3. plot of observed vs fitted values
    
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predictsSites$Species_richness,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    
    ## 4. Check for spatial autocorrelation
    
    sa_test<-SpatialAutocorrelationTest(model=get(x), ranefGrouping = "SS")
    
    #summary(sa_test)
    
    # percentage of studies that show spatial autocorrelation?
    perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100
    
    
    sa_test_vals <- as.data.frame(sa_test$P)
    sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)
    
    label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")
    
    p4 <- ggplot(data = sa_test_vals ) +
      geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
      geom_vline(xintercept = 0.05, col = "red") +
      geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
      theme_bw() +
      ylim(c(0, 100)) +
      xlab("P-value") +
      ylab("Frequency") +
      theme(panel.grid = element_blank(), 
            aspect.ratio = 1)
    
    cowplot::plot_grid(p2, p3, p4,
                       labels = c("(a)", "(b)", "(c)"))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p2, p3, p4)
    
  }
  
}



