
############################################################
#                                                          #
#    Additional tests: Chao-estimated species richness     #
#                                                          #
############################################################

# in this script, models are run using Chao-estimated species richness

#### Set up ####
rm(list = ls())

# set up directories
# dataDir <- "Data/"
inDir<- "5_RunLUIClimateModels/"
outDir <- "10_Additional_Tests/Chao_SpeciesRichness/"
if(!dir.exists(outDir)) dir.create(outDir)

# Load required libraries
library(StatisticalModels)
library(dplyr)
library(ggplot2)
library(cowplot)

# source in additional functions
source("0_Functions.R")

# load in the dataset 
sites <- readRDS(file = paste0(inDir, "PREDICTSSitesClimate_Data.rds"))

# round the estimated species richness values to integers.
sites$ChaoR <- round(sites$ChaoR,0)

# remove sites that do not have an estimate of ChaoR
sites <- sites[!is.na(sites$ChaoR), ]
# 5152 rows

# assess coverage of the data
table(sites$Order, sites$LUI) 

#             Primary vegetation Secondary vegetation Agriculture_Low Agriculture_High
# Coleoptera                 446                  232             169              282
# Diptera                     76                   95              54              111
# Hemiptera                   77                   83             196               99
# Hymenoptera                367                  410             441             1125
# Lepidoptera                217                  384             207               81


length(unique(sites$SSBS)) # 4058 sites
length(unique(sites$SS)) # 154 studies

############################################################
#                                                          #
####      Run the model with ChaoR as the response      ####
#                                                          #
############################################################

# 1. Chao Richness, mean anomaly)

MeanAnomalyModelChaoR <- GLMER(modelData = sites,
                               responseVar = "ChaoR",
                               fitFamily = "poisson",
                               fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("Total_abundance","SSBS"), 
                               maxIters = 20000)


save(MeanAnomalyModelChaoR, file = paste0(outDir, "/MeanAnomalyModelChaoR.rdata"))

summary(MeanAnomalyModelChaoR$model)



############################################################
#                                                          #
####                 Plot the results                   ####
#                                                          #
############################################################


#### 1. mean anom ChaoR ####

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS),
                        length.out = 150),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelChaoR$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = sites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$ChaoR <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

# coleoptera
C_QPV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Primary vegetation" & MeanAnomalyModelChaoR$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QSV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Secondary vegetation" & MeanAnomalyModelChaoR$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAL <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_Low" & MeanAnomalyModelChaoR$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAH <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_High" & MeanAnomalyModelChaoR$data$Order == "Coleoptera"],
  probs = exclQuantiles)

# diptera
D_QPV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Primary vegetation" & MeanAnomalyModelChaoR$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QSV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Secondary vegetation" & MeanAnomalyModelChaoR$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAL <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_Low" & MeanAnomalyModelChaoR$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAH <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_High" & MeanAnomalyModelChaoR$data$Order == "Diptera"],
  probs = exclQuantiles)

# hemiptera
He_QPV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Primary vegetation" & MeanAnomalyModelChaoR$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QSV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Secondary vegetation" & MeanAnomalyModelChaoR$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAL <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_Low" & MeanAnomalyModelChaoR$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAH <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_High" & MeanAnomalyModelChaoR$data$Order == "Hemiptera"],
  probs = exclQuantiles)

# hymenoptera
Hy_QPV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Primary vegetation" & MeanAnomalyModelChaoR$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QSV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Secondary vegetation" & MeanAnomalyModelChaoR$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAL <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_Low" & MeanAnomalyModelChaoR$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAH <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_High" & MeanAnomalyModelChaoR$data$Order == "Hymenoptera"],
  probs = exclQuantiles)

# lepidoptera
L_QPV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Primary vegetation" & MeanAnomalyModelChaoR$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QSV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Secondary vegetation" & MeanAnomalyModelChaoR$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAL <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_Low" & MeanAnomalyModelChaoR$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAH <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$LUI=="Agriculture_High" & MeanAnomalyModelChaoR$data$Order == "Lepidoptera"],
  probs = exclQuantiles)


# predict the results
c.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelChaoR$model, data = nd)

# back transform the abundance values
c.preds.tmean <- exp(c.preds.tmean)

# create list of matrices
number_of_chunks = 5
list_c.preds.tmean <- lapply(seq(1, NROW(c.preds.tmean), ceiling(NROW(c.preds.tmean)/number_of_chunks)),
                             function(i) c.preds.tmean[i:min(i + ceiling(NROW(c.preds.tmean)/number_of_chunks) - 1, NROW(c.preds.tmean)),])

# name them
names(list_c.preds.tmean) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")

# sweep according to the refRow
list_c.preds.tmean <- lapply(list_c.preds.tmean,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[refRow[1],],FUN="/") 
})

list2env(list_c.preds.tmean,globalenv())

# split nd by order
Order<- paste0("nd_",nd$Order)
# create a list of data frames
by_Order <- split(nd,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QPV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QPV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QSV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QSV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QAL[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QAL[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QAH[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QAH[2])),] <- NA

Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > D_QPV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < D_QPV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < D_QSV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > D_QSV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS < D_QAL[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS > D_QAL[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS < D_QAH[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS > D_QAH[2])),] <- NA

Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QPV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QPV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QSV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QSV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QAL[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QAL[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QAH[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QAH[2])),] <- NA

Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QPV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QPV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QSV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QSV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QAL[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QAL[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QAH[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QAH[2])),] <- NA

Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QPV[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QPV[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Secondary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QSV[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Secondary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QSV[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_Low") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QAL[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_Low") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QAL[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_High") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QAH[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_High") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = median,na.rm=TRUE))*100)-100
nd_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = median,na.rm=TRUE))*100)-100
nd_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                  FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                  FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100



# set factor levels
# nd$LUI <- factor(nd$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Coleoptera$LUI <- factor(nd_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Diptera$LUI <- factor(nd_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hemiptera$LUI <- factor(nd_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hymenoptera$LUI <- factor(nd_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Lepidoptera$LUI <- factor(nd_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Coleoptera$PredLower, ymax = nd_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500), limits = c(-100, 550)) +
  ylab("Change in Chao species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(linewidth = 0.2)) + 
  ggtitle("Coleoptera")

p_diptera <- ggplot(data = nd_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Diptera$PredLower, ymax = nd_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500), limits = c(-100, 550)) +
  ylab("Change in Chao species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(linewidth = 0.2)) + 
  ggtitle("Diptera")

p_hemiptera <- ggplot(data = nd_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Hemiptera$PredLower, ymax = nd_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500), limits = c(-100, 550)) +
  ylab("Change in Chao species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(linewidth = 0.2)) + 
  ggtitle("Hemiptera")

p_hymenoptera <- ggplot(data = nd_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Hymenoptera$PredLower, ymax = nd_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500), limits = c(-100, 550)) +
  ylab("Change in Chao species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(linewidth = 0.2)) + 
  ggtitle("Hymenoptera")

p_lepidoptera <- ggplot(data = nd_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Lepidoptera$PredLower, ymax = nd_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500), limits = c(-100, 550)) +
  ylab("Change in Chao species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(linewidth = 0.2)) + 
  ggtitle("Lepidoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 11), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomChao <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera)
MeanAnomChao <- cowplot::plot_grid(MeanAnomChao,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
ggsave(filename = paste0(outDir, "MeanAnomChao.pdf"), plot = MeanAnomChao, width = 200, height = 150, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "MeanAnomChao.jpeg"), plot = MeanAnomChao, width = 200, height = 150, units = "mm", dpi = 300)

