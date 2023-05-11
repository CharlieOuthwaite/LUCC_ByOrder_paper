##%######################################################%##
#                                                          #
####             Separate models per order              ####
#                                                          #
##%######################################################%##

# rather than including order as an interactive term, we run a
# separate model for each insect order using a subset of the data. 

# clear working environment
rm(list = ls())

# set directories 
datadir <- "5_RunLUIClimateModels/"
outdir <- "10_Additional_Tests/Order_level_models/"
if(!dir.exists(outdir)) dir.create(outdir)

source("0_Functions.R")


# load in the data
predictsSites <- readRDS(file = paste0(datadir,"PREDICTSSitesClimate_Data.rds"))


#### loop through and run a model for each order ####

orders <- c("Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera")

# order <- orders[1]

for(order in orders){
  
  order_data <- predictsSites[predictsSites$Order == order, ]
  
  
  model_data <- order_data[!is.na(order_data$LogAbund), ]
  model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]
  

  # run  model
  MeanAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                 fixedStruct = "LUI * StdTmeanAnomalyRS",
                                 randomStruct = "(1|SS)+(1|SSB)",
                                 saveVars = c("SSBS"))
  
  # save the model output
  save(MeanAnomalyModelAbund, file = paste0(outdir, "MeanAnomalyModelAbund_", order, ".rdata"))
  
  
  # ii. Richness, mean anomaly
  
  model_data <- order_data[!is.na(order_data$StdTmeanAnomalyRS), ]
  
  MeanAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                fixedStruct = "LUI * StdTmeanAnomalyRS",
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                saveVars = c("SSBS"))
  
  summary(MeanAnomalyModelRich$model)
  
  save(MeanAnomalyModelRich, file = paste0(outdir, "MeanAnomalyModelRich_", order, ".rdata"))
  
  
  
}


#### Loop through orders and plot results ####


for(order in orders){
  
  # load the models
  
  load(file = paste0(outdir, "MeanAnomalyModelAbund_", order, ".rdata"))
  load(file = paste0(outdir, "MeanAnomalyModelRich_", order, ".rdata"))
  
  
  # set quantiles of predicted result to be presented in the plots
  exclQuantiles <- c(0.025,0.975)
  

  nd <- expand.grid(
    StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                          to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                          length.out = 100),
    LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
               levels = levels(MeanAnomalyModelAbund$data$LUI)))
  
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
  
  QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
    MeanAnomalyModelAbund$data$LUI=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
    MeanAnomalyModelAbund$data$LUI=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
    MeanAnomalyModelAbund$data$LUI=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
    MeanAnomalyModelAbund$data$LUI=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)
  
  # back transform the abundance values
  a.preds.tmean <- exp(a.preds.tmean)-0.01
  
  # convert to relative to reference
  a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')
  
  # remove anything above and below the quantiles
  a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
  a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
  a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
  a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
  a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
  a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
  a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
  a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA
  
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
    ggtitle(order)
  
  
  ## now the species richness plot 
  nd2 <- expand.grid(
    StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                          to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                          length.out = 100),
    LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
               levels = levels(MeanAnomalyModelRich$data$LUI)))
  
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
  
  QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
    MeanAnomalyModelRich$data$LUI=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
    MeanAnomalyModelRich$data$LUI=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
    MeanAnomalyModelRich$data$LUI=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
    MeanAnomalyModelRich$data$LUI=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd2, nIters = 10000)
  
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
    scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100), limits = c(-100, 100)) +
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
    ggtitle(order)
  
  
  # combine plots
  cowplot::plot_grid(p1, p2)
  
  ggsave(filename = paste0(outdir, order, "_Abun_Rich.png"), plot = last_plot(), width = 183, height = 100, units = "mm", dpi = 300)
  
  
  
  
}

# check R-squared values

R2GLMER (MeanAnomalyModelAbund_Coleoptera$model)
# $conditional
# [1] 0.5924863
# 
# $marginal
# [1] 0.04200098

R2GLMER (MeanAnomalyModelAbund_Diptera$model)
# $conditional
# [1] 0.5303029
# 
# $marginal
# [1] 0.06720692

R2GLMER (MeanAnomalyModelAbund_Hemiptera$model)
# $conditional
# [1] 0.4286731
# 
# $marginal
# [1] 0.07476181

R2GLMER (MeanAnomalyModelAbund_Hymenoptera$model)
# $conditional
# [1] 0.4851417
# 
# $marginal
# [1] 0.01828185

R2GLMER (MeanAnomalyModelAbund_Lepidoptera$model)
# $conditional
# [1] 0.5404108
# 
# $marginal
# [1] 0.03082085

R2GLMER (MeanAnomalyModelRich_Coleoptera$model)
# $conditional
# [1] 0.7595814
# 
# $marginal
# [1] 0.01589639

R2GLMER (MeanAnomalyModelRich_Diptera$model)
# $conditional
# [1] 0.8044178
# 
# $marginal
# [1] 0.1016308

R2GLMER (MeanAnomalyModelRich_Hemiptera$model)
# $conditional
# [1] 0.6736128
# 
# $marginal
# [1] 0.06802885

R2GLMER (MeanAnomalyModelRich_Hymenoptera$model)
# $conditional
# [1] 0.6459363
# 
# $marginal
# [1] 0.01125481

R2GLMER (MeanAnomalyModelRich_Lepidoptera$model)
# $conditional
# [1] 0.8252371
# 
# $marginal
# [1] 0.007577729