
############################################################
#                                                          #
#    Additional tests: Chao-estimated species richness     #
#                                                          #
############################################################

# in this script, models are run using Chao-estimated species richness

#### Set up ####

# set up directories
dataDir <- "Data/"
inDir<- "4_PREDICTSMatchClimateIndex/"
outDir <- "10_Additional_Tests/Chao_SpeciesRichness/"
if(!dir.exists(outDir)) dir.create(outDir)

# Load required libraries
library(predictsFunctions)
library(dplyr)
library(ggplot2)
library(cowplot)

# source in additional functions
source("0_Functions.R")

# Set the path to copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)
# 3250404 obs. of 67 variables

# Select only data for insects
predicts <- predicts[(predicts$Class=="Insecta"),]
# 935078 obs. of 67 variables

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
# Correcting 0 missing sampling effort values
# Re-scaling sampling effort
# Correcting 870378 values for sensitivity to sampling effort # matches Outhwaite et al.

# check diversity metrics
table(predicts$Diversity_metric)

# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those entries that are the problem
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]
predicts <- droplevels(predicts)
# 934845 obs. of 67 variables

# MergeSites

# Merge sites that have the same coordinates, from the same study and same taxonomic family (e.g. multiple traps on a single transect)
predicts <- MergeSites(diversity = predicts)

# prepare PREDICTS
# remove entries without Order
predicts <- droplevels(predicts[(predicts$Order!=""),])
# 826016 obs. of 67 variables

# filter only for the orders with at least 600 unique sites
predicts <- predicts %>% filter(Order %in% c("Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera", "Hemiptera")) %>% droplevels()
# 814738 obs. of 67 variables

# convert Order to a "factor"
predicts$Order <- as.factor(predicts$Order)

# check
table(predicts$Order)

# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera   
#     390020       27157       49119      176639      167464      

# Split predicts into separate data frames according to insect Order
OrderName <- paste0("",predicts$Order)

by_Order <- split(predicts,OrderName)

# extract data frames from list into global environment
list2env(by_Order,globalenv())

# Calculate site metrics of diversity for each order, include extra columns:
# Predominant_land_use
# SSB
# SSBS
# Biome
# Order
# species richness estimators

# droplevels() drops unused factor levels. 
# This is particularly useful if we want to drop factor levels that are no longer 
# used due to subsetting a vector or a data frame (as we did with split()). 

Coleoptera <- droplevels(Coleoptera)
Coleoptera <- SiteMetrics(diversity = Coleoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order"),
                          srEstimators = "Chao")
Diptera <- droplevels(Diptera)
Diptera <- SiteMetrics(diversity = Diptera,
                       extra.cols = c("Predominant_land_use",
                                      "SSB","SSBS", "Biome","Order"),
                       srEstimators = "Chao")
Hemiptera <- droplevels(Hemiptera)
Hemiptera <- SiteMetrics(diversity = Hemiptera,
                         extra.cols = c("Predominant_land_use",
                                        "SSB","SSBS", "Biome","Order"),
                         srEstimators = "Chao")
Hymenoptera <- droplevels(Hymenoptera)
Hymenoptera <- SiteMetrics(diversity = Hymenoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order"),
                           srEstimators = "Chao")
Lepidoptera <- droplevels(Lepidoptera)
Lepidoptera <- SiteMetrics(diversity = Lepidoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order"),
                           srEstimators = "Chao")


# merge all sites_Order data frames into one called "sites" using rbind()
sites <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera)
# 10746 obs. of 23 variables

# Chao cannot be estimated for all sites, so drop those that it can't and see what is left
sites <- sites[!is.na(sites$ChaoR),] 
# 6585 obs. of 23 variables

# simple plot of estimated species richness against sampled species richness
jpeg(file="10_Additional_Tests/Species Richness vs Chao.jpeg")

par(tck=-0.01,mgp=c(1.6,0.2,0),mar=c(2.7,2.7,0.2,0.2),las=1)
plot(sites$Species_richness+1,sites$ChaoR+1,log="xy",pch=16,
     xlab="Sampled species richness",ylab="Estimated species richness")
abline(0,1,lwd=2,col="#ff0000")

dev.off()

# histogram of the ratios of estimated to sampled species richness:
jpeg(file="10_Additional_Tests/Estimated to Sampled Species Richness Ratio.jpeg")

par(tck=-0.01,mgp=c(1.6,0.2,0),mar=c(2.7,2.7,0.2,0.2),las=1)
hist(x = log10((sites$ChaoR-sites$Species_richness)+1),
     xaxt="n",xlab="Estimated - Sampled richness",main=NULL)
axis(1,at=log10(c(0,1,5,10,100,1000)+1),labels=c(0,1,5,10,100,1000))

dev.off()

# Test whether completeness of sampling is related to sampled or estimated species richness:
jpeg(file="10_Additional_Tests/Completeness of Sampling.jpeg")

par(mfrow=c(1,2),tck=-0.01,mgp=c(1.6,0.2,0),mar=c(2.7,2.7,0.2,0.2),las=1)
plot(x = sites$Species_richness+1,y = (sites$ChaoR - sites$Species_richness)+1,log="xy",pch=16,
     xlab="Sampled species richness",ylab="Estimated - Sampled richness",
     xaxt="n",yaxt="n")
axis(1,at=c(1,11,51,101,501),labels=c(0,10,50,100,500))
axis(2,at=c(1,11,51,101,501),labels=c(0,10,50,100,500))
plot(x = sites$ChaoR+1,y = (sites$ChaoR - sites$Species_richness)+1,log="xy",pch=16,
     xlab="Estimated species richness",ylab="Estimated - Sampled richness",
     xaxt="n",yaxt="n")
axis(1,at=c(1,11,51,101,501,1001),labels=c(0,10,50,100,500,1000))
axis(2,at=c(1,11,51,101,501),labels=c(0,10,50,100,500))

dev.off()

# round the estimated species richness values to integers.
sites$ChaoR <- round(sites$ChaoR,0)

# First, we will rearrange the land-use classification a bit
# rename "Predominant_land_use" to "LandUse"
sites$LandUse <- paste(sites$Predominant_land_use)

# Drop classification where land use could not be identified
sites$LandUse[(sites$LandUse=="Cannot decide")] <- NA

# drop classification where use intensity could not be identified
sites$Use_intensity[sites$Use_intensity=="Cannot decide"] <- NA

# Now make the variable a factor, and set the reference level to primary vegetation
sites$LandUse <- factor(sites$LandUse)
sites$LandUse <- relevel(sites$LandUse,ref="Primary vegetation")

# combine LandUse (LU) and Use Intensity (UI) into new variable Land Use Intensity (LUI)
sites$LUI <- paste0(sites$LandUse,'_',sites$Use_intensity)
sites$LUI[grep("NA",sites$LUI)] <- NA # where "NA" appears in the UI field, drop classification

# recode according to land use and use intensity combinations
sites$LUI <- dplyr::recode(sites$LUI,
                           'Primary vegetation_Minimal use' = 'Primary vegetation',
                           'Cropland_Light use' = 'Agriculture_High',
                           'Secondary vegetation (indeterminate age)_Minimal use' = 'Secondary vegetation',
                           'Urban_Light use' = 'Urban',
                           'Secondary vegetation (indeterminate age)_Light use' = 'Secondary vegetation',
                           'Cropland_Intense use' = 'Agriculture_High',
                           'Cropland_Minimal use' = 'Agriculture_Low',
                           'Pasture_Light use' = 'Agriculture_Low',
                           'Pasture_Minimal use' = 'Agriculture_Low',
                           'Intermediate secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Secondary vegetation (indeterminate age)_Intense use' = 'Secondary vegetation',
                           'Pasture_Intense use' = 'Agriculture_High',
                           'Urban_Minimal use' = 'Urban',
                           'Primary vegetation_Light use' = 'Primary vegetation',
                           'Young secondary vegetation_Light use' = 'Secondary vegetation',
                           'Urban_Intense use' = 'Urban',
                           'Primary vegetation_Intense use' = 'Primary vegetation',
                           'Young secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Minimal use' = 'Agriculture_Low',
                           'Plantation forest_Intense use' = 'Agriculture_High',
                           'Young secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Light use' = 'Agriculture_High',
                           'Mature secondary vegetation_Light use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Light use' = 'Secondary vegetation')

# 6948 obs. of 25 variables

sites$Use_intensity[((sites$LandUse=="Mature secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Intermediate secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Young secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"

# remove the urban sites and sites that are NA in LUI
sites <- sites[!sites$LUI == "Urban", ]
sites <- sites[!is.na(sites$LUI), ]

# 5920 obs. of 25 variables
sites <- droplevels(sites)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]
# 5920 obs. of 25 variables

# load in original dataset and match climate anomaly data to sites
predictsSites <- readRDS(paste0(inDir,"PREDICTSSites_Climate.rds"))
predictsSites <- predictsSites@data
# 8884 obs. of 35 variables

# keep only SSBS and climate anomaly variables
predictsSites <- predictsSites[ , c("SSBS", "StdTmeanAnomaly", "StdTmaxAnomaly")]
# 8884 obs. of 3 variables

# Select unique sites only so there is no duplication
predictsSites <- unique(predictsSites)

# combine to get the anomaly info
sites2 <- merge(sites, predictsSites, by = "SSBS")
# 5920 obs. of 27 variables

sites2$LUI <- factor(sites2$LUI)
sites2$LUI <- relevel(sites2$LUI,ref="Primary vegetation")

# organise the climate anomaly data
sites2$StdTmeanAnomalyRS <- StdCenterPredictor(sites2$StdTmeanAnomaly)

# rescale the variable
sites2$StdTmaxAnomalyRS <- StdCenterPredictor(sites2$StdTmaxAnomaly)


# drop levels
sites2 <- droplevels(sites2)

save(sites2, file = paste0(outDir, "PREDICTS_sites_ChaoR.rdata"))


table(sites2$Order, sites2$LUI)

#             Primary vegetation Agriculture_High Agriculture_Low Secondary vegetation
# Coleoptera                 512              282             228                  339
# Diptera                     77              143              86                   95
# Hemiptera                  143               99             255                  110
# Hymenoptera                367             1162             478                  423
# Lepidoptera                283               81             266                  491

############################################################
#                                                          #
####      Run the model with ChaoR as the response      ####
#                                                          #
############################################################


# 1. Chao Richness, mean anomaly)

MeanAnomalyModelChaoR <- GLMER(modelData = sites2,responseVar = "ChaoR",
                              fitFamily = "poisson",
                              fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("Total_abundance","SSBS"))

save(MeanAnomalyModelChaoR, file = paste0(outDir, "/MeanAnomalyModelChaoR.rdata"))


summary(MeanAnomalyModelChaoR$model)
# ChaoR ~ LUI * StdTmeanAnomalyRS * Order + (1 | SS) + (1 | SSB) + (1 | SSBS)


# # 1. Chao Richness, max anomaly
# MaxAnomalyModelChaoR <- GLMER(modelData = sites2,responseVar = "ChaoR",
#                               fitFamily = "poisson",
#                               fixedStruct = "LUI * StdTmaxAnomalyRS * Order",
#                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
#                               saveVars = c("Total_abundance","SSBS"))
# 
# 
# save(MaxAnomalyModelChaoR, file = paste0(outDir, "MaxAnomalyModelChaoR.rdata"))
# 
# summary(MaxAnomalyModelChaoR$model)
# # ChaoR ~ LUI + poly(StdTmaxAnomalyRS, 1) + LUI:poly(StdTmaxAnomalyRS,      1) + (1 | SS) + (1 | SSB) + (1 | SSBS)


############################################################
#                                                          #
####                 Plot the results                   ####
#                                                          #
############################################################


# load models
load(paste0(outDir, "MeanAnomalyModelChaoR.rdata"))
# load(paste0(outDir, "MaxAnomalyModelChaoR.rdata"))

# load sites2
load(paste0(outDir, "PREDICTS_sites_ChaoR.rdata"))

#### 1. mean anom ChaoR ####

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelChaoR$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = sites2$StdTmeanAnomaly)

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
c.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelChaoR$model,data = nd)

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
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500, 600, 700, 800, 900), limits = c(-100, 950)) +
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
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Coleoptera")

p_diptera <- ggplot(data = nd_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Diptera$PredLower, ymax = nd_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500, 600), limits = c(-100, 600)) +
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
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Diptera")

p_hemiptera <- ggplot(data = nd_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Hemiptera$PredLower, ymax = nd_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500, 600), limits = c(-100, 600)) +
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
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hemiptera")

p_hymenoptera <- ggplot(data = nd_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Hymenoptera$PredLower, ymax = nd_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500, 600), limits = c(-100, 600)) +
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
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hymenoptera")

p_lepidoptera <- ggplot(data = nd_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), linewidth = 0.75) +
  geom_ribbon(aes(ymin = nd_Lepidoptera$PredLower, ymax = nd_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100,  0, 100, 200, 300, 400, 500, 600), limits = c(-100, 600)) +
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
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
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

# #### maximum anomaly ####
# 
# nd2 <- expand.grid(
#   StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS),
#                         to = max(MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS),
#                         length.out = 100),
#   LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(MaxAnomalyModelChaoR$data$LUI)),
#   Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))
# 
# # back transform the predictor
# nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd2$StdTmaxAnomalyRS,
#   originalX = sites2$StdTmaxAnomaly)
# 
# # set richness and abundance to 0 - to be predicted
# nd2$ChaoR <- 0
# 
# # reference for % difference = primary vegetation and positive anomaly closest to 0
# refRow <- which((nd2$LUI=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))))
# 
# # adjust plot 1: mean anomaly and abundance
# 
# # coleoptera
# C_QPV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Primary vegetation" & MaxAnomalyModelChao$data$Order == "Coleoptera"],
#   probs = exclQuantiles)
# C_QSV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Secondary vegetation" & MaxAnomalyModelChao$data$Order == "Coleoptera"],
#   probs = exclQuantiles)
# C_QAL <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_Low" & MaxAnomalyModelChao$data$Order == "Coleoptera"],
#   probs = exclQuantiles)
# C_QAH <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_High" & MaxAnomalyModelChao$data$Order == "Coleoptera"],
#   probs = exclQuantiles)
# 
# # diptera
# D_QPV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Primary vegetation" & MaxAnomalyModelChao$data$Order == "Diptera"],
#   probs = exclQuantiles)
# D_QSV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Secondary vegetation" & MaxAnomalyModelChao$data$Order == "Diptera"],
#   probs = exclQuantiles)
# D_QAL <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_Low" & MaxAnomalyModelChao$data$Order == "Diptera"],
#   probs = exclQuantiles)
# D_QAH <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_High" & MaxAnomalyModelChao$data$Order == "Diptera"],
#   probs = exclQuantiles)
# 
# # hemiptera
# He_QPV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Primary vegetation" & MaxAnomalyModelChao$data$Order == "Hemiptera"],
#   probs = exclQuantiles)
# He_QSV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Secondary vegetation" & MaxAnomalyModelChao$data$Order == "Hemiptera"],
#   probs = exclQuantiles)
# He_QAL <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_Low" & MaxAnomalyModelChao$data$Order == "Hemiptera"],
#   probs = exclQuantiles)
# He_QAH <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_High" & MaxAnomalyModelChao$data$Order == "Hemiptera"],
#   probs = exclQuantiles)
# 
# # hymenoptera
# Hy_QPV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Primary vegetation" & MaxAnomalyModelChao$data$Order == "Hymenoptera"],
#   probs = exclQuantiles)
# Hy_QSV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Secondary vegetation" & MaxAnomalyModelChao$data$Order == "Hymenoptera"],
#   probs = exclQuantiles)
# Hy_QAL <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_Low" & MaxAnomalyModelChao$data$Order == "Hymenoptera"],
#   probs = exclQuantiles)
# Hy_QAH <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_High" & MaxAnomalyModelChao$data$Order == "Hymenoptera"],
#   probs = exclQuantiles)
# 
# # lepidoptera
# L_QPV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Primary vegetation" & MaxAnomalyModelChao$data$Order == "Lepidoptera"],
#   probs = exclQuantiles)
# L_QSV <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Secondary vegetation" & MaxAnomalyModelChao$data$Order == "Lepidoptera"],
#   probs = exclQuantiles)
# L_QAL <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_Low" & MaxAnomalyModelChao$data$Order == "Lepidoptera"],
#   probs = exclQuantiles)
# L_QAH <- quantile(x = MaxAnomalyModelChao$data$StdTmeanAnomalyRS[
#   MaxAnomalyModelChao$data$LUI=="Agriculture_High" & MaxAnomalyModelChao$data$Order == "Lepidoptera"],
#   probs = exclQuantiles)
# 
# # predict the results
# c.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelChaoR$model,data = nd2)
# 
# # back transform the abundance values
# c.preds.tmax <- exp(c.preds.tmax)
# 
# # create list of matrices
# number_of_chunks = 5
# list_c.preds.tmax <- lapply(seq(1, NROW(c.preds.tmax), ceiling(NROW(c.preds.tmax)/number_of_chunks)),
#                              function(i) c.preds.tmax[i:min(i + ceiling(NROW(c.preds.tmax)/number_of_chunks) - 1, NROW(c.preds.tmax)),])
# 
# # name them
# names(list_c.preds.tmax) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")
# 
# # sweep according to the refRow 
# list_c.preds.tmax <- lapply(list_c.preds.tmax,FUN=function(x){
#   sweep (x=x, MARGIN = 2, STATS=x[refRow[1],],FUN="/") 
# })
# 
# list2env(list_c.preds.tmax,globalenv())
# 
# # split nd by order
# Order<- paste0("nd2_",nd2$Order)
# # create a list of data frames
# by_Order <- split(nd2,Order)
# list2env(by_Order,globalenv())
# 
# # remove anything above and below the quantiles
# 
# Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QPV[2])),] <- NA
# Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QPV[1])),] <- NA
# Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QSV[1])),] <- NA
# Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QSV[2])),] <- NA
# Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QAL[1])),] <- NA
# Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QAL[2])),] <- NA
# Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QAH[1])),] <- NA
# Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QAH[2])),] <- NA
# 
# Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > D_QPV[2])),] <- NA
# Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < D_QPV[1])),] <- NA
# Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < D_QSV[1])),] <- NA
# Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > D_QSV[2])),] <- NA
# Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS < D_QAL[1])),] <- NA
# Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS > D_QAL[2])),] <- NA
# Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS < D_QAH[1])),] <- NA
# Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS > D_QAH[2])),] <- NA
# 
# Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QPV[2])),] <- NA
# Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QPV[1])),] <- NA
# Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QSV[1])),] <- NA
# Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QSV[2])),] <- NA
# Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QAL[1])),] <- NA
# Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QAL[2])),] <- NA
# Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QAH[1])),] <- NA
# Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QAH[2])),] <- NA
# 
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QPV[2])),] <- NA
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QPV[1])),] <- NA
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QSV[1])),] <- NA
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QSV[2])),] <- NA
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QAL[1])),] <- NA
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QAL[2])),] <- NA
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QAH[1])),] <- NA
# Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QAH[2])),] <- NA
# 
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QPV[2])),] <- NA
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QPV[1])),] <- NA
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QSV[1])),] <- NA
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QSV[2])),] <- NA
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QAL[1])),] <- NA
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QAL[2])),] <- NA
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QAH[1])),] <- NA
# Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QAH[2])),] <- NA
# 
# # Get the median, upper and lower quants for the plot
# 
# nd2_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
#                                     FUN = median,na.rm=TRUE))*100)-100
# nd2_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
#                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
#                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd2_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
#                                  FUN = median,na.rm=TRUE))*100)-100
# nd2_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
#                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
#                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd2_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
#                                    FUN = median,na.rm=TRUE))*100)-100
# nd2_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
#                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
#                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd2_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
#                                      FUN = median,na.rm=TRUE))*100)-100
# nd2_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
#                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
#                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# nd2_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
#                                      FUN = median,na.rm=TRUE))*100)-100
# nd2_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
#                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
#                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# # set factor levels
# nd2_Coleoptera$LUI <- factor(nd2_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# nd2_Diptera$LUI <- factor(nd2_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# nd2_Hemiptera$LUI <- factor(nd2_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# nd2_Hymenoptera$LUI <- factor(nd2_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# nd2_Lepidoptera$LUI <- factor(nd2_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
# 
# # plot
# 
# p_coleoptera <- ggplot(data = nd2_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = LUI), size = 0.75) +
#   geom_ribbon(aes(ymin = nd2_Coleoptera$PredLower, ymax = nd2_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
#   geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   theme_bw() + 
#   scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
#   #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
#   scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
#   ylab("Change in Chao species richness (%)") +
#   xlab("Standardised Maximum Temperature Anomaly") +
#   #xlim(c(-1, 5)) +
#   #ylim(c(-65, 60)) + 
#   theme(aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 7),
#         legend.position = "none",
#         #legend.position = c(0.2, 0.8),
#         #legend.background = element_blank(), 
#         #legend.text = element_text(size = 6), 
#         #legend.title = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(size = 0.2),
#         panel.border = element_rect(size = 0.2), 
#         axis.ticks = element_line(size = 0.2)) + 
#   ggtitle("Coleoptera")
# 
# p_diptera <- ggplot(data = nd2_Diptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = LUI), size = 0.75) +
#   geom_ribbon(aes(ymin = nd2_Diptera$PredLower, ymax = nd2_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
#   geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   theme_bw() + 
#   scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
#   #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
#   scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
#   ylab("Change in Chao species richness (%)") +
#   xlab("Standardised Maximum Temperature Anomaly") +
#   #xlim(c(-1, 5)) +
#   #ylim(c(-65, 60)) + 
#   theme(aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 7),
#         legend.position = "none",
#         #legend.position = c(0.2, 0.8),
#         #legend.background = element_blank(), 
#         #legend.text = element_text(size = 6), 
#         #legend.title = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(size = 0.2),
#         panel.border = element_rect(size = 0.2), 
#         axis.ticks = element_line(size = 0.2)) + 
#   ggtitle("Diptera")
# 
# p_hemiptera <- ggplot(data = nd2_Hemiptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = LUI), size = 0.75) +
#   geom_ribbon(aes(ymin = nd2_Hemiptera$PredLower, ymax = nd2_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
#   geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   theme_bw() + 
#   scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
#   #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
#   scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
#   ylab("Change in Chao species richness (%)") +
#   xlab("Standardised Maximum Temperature Anomaly") +
#   #xlim(c(-1, 5)) +
#   #ylim(c(-65, 60)) + 
#   theme(aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 7),
#         legend.position = "none",
#         #legend.position = c(0.2, 0.8),
#         #legend.background = element_blank(), 
#         #legend.text = element_text(size = 6), 
#         #legend.title = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(size = 0.2),
#         panel.border = element_rect(size = 0.2), 
#         axis.ticks = element_line(size = 0.2)) + 
#   ggtitle("Hemiptera")
# 
# p_hymenoptera <- ggplot(data = nd2_Hymenoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = LUI), size = 0.75) +
#   geom_ribbon(aes(ymin = nd2_Hymenoptera$PredLower, ymax = nd2_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
#   geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   theme_bw() + 
#   scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
#   #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
#   scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
#   ylab("Change in Chao species richness (%)") +
#   xlab("Standardised Maximum Temperature Anomaly") +
#   #xlim(c(-1, 5)) +
#   #ylim(c(-65, 60)) + 
#   theme(aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 7),
#         legend.position = "none",
#         #legend.position = c(0.2, 0.8),
#         #legend.background = element_blank(), 
#         #legend.text = element_text(size = 6), 
#         #legend.title = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(size = 0.2),
#         panel.border = element_rect(size = 0.2), 
#         axis.ticks = element_line(size = 0.2)) + 
#   ggtitle("Hymenoptera")
# 
# p_lepidoptera <- ggplot(data = nd2_Lepidoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = LUI), size = 0.75) +
#   geom_ribbon(aes(ymin = nd2_Lepidoptera$PredLower, ymax = nd2_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
#   geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
#   scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
#   theme_bw() + 
#   scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
#   #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
#   scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
#   ylab("Change in Chao species richness (%)") +
#   xlab("Standardised Maximum Temperature Anomaly") +
#   #xlim(c(-1, 5)) +
#   #ylim(c(-65, 60)) + 
#   theme(aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 7),
#         legend.position = "none",
#         #legend.position = c(0.2, 0.8),
#         #legend.background = element_blank(), 
#         #legend.text = element_text(size = 6), 
#         #legend.title = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(size = 0.2),
#         panel.border = element_rect(size = 0.2), 
#         axis.ticks = element_line(size = 0.2)) + 
#   ggtitle("Lepidoptera")
# 
# # get the legend
# legend <- get_legend(
#   p_coleoptera +
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom",
#           legend.background = element_blank(), 
#           legend.text = element_text(size = 11), 
#           legend.title = element_blank())
# )
# 
# 
# # put them all together to save them
# MaxAnomChao <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera)
# MaxAnomChao <- cowplot::plot_grid(MaxAnomChao,legend,ncol=1, rel_heights = c(1,0.1))
# 
# # save the ggplots
# ggsave(filename = paste0(outDir, "MaxAnomChao.pdf"), plot = MaxAnomChao, width = 200, height = 150, units = "mm", dpi = 300)
