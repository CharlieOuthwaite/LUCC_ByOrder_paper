############################################################
#                                                          #
#           Jack-knifing removal of each family            #
#                                                          #
############################################################

# Here we remove each family from the dataset before rerunning the 
# analysis to see if any one family has a strong influence on the results. 
# So only one family is removed each time, and then added back in for 
# subsequent tests. 

rm(list = ls())

# libraries
library(StatisticalModels)
library(ggplot2)
library(cowplot)
library(dplyr)

# directories
dataDir <- "Data/"
outdir <- "10_Additional_Tests/Jackknife_Families/"
if(!dir.exists(outdir)) dir.create(outdir)

# source in additional functions
source("0_Functions.R")


#### 1. Load and check dataset ####

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

# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those entries that are the problem
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]
# 934845 obs. of 67 variables

# Merge sites that have the same coordinates, from the same study and same taxonomic family (e.g. multiple traps on a single transect)
predicts <- MergeSites(diversity = predicts)
# 826292 obs. of 67 variables

# remove entries without Order
predicts <- droplevels(predicts[(predicts$Order!=""),])
# 826016 obs. of 67 variables

# convert Order to a "factor"
predicts$Order <- as.factor(predicts$Order)

# remove orders not included
predicts <- predicts[predicts$Order %in% c("Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"), ]


# how many unique families?
fams <- droplevels(unique(predicts$Family))
length(fams) # 357

predicts <- droplevels(predicts)

# how many famliies per order
results <- predicts %>%
  group_by(Order) %>%
  summarize(count = n_distinct(Family))

# # A tibble: 5 × 2
# Order       count
# <fct>       <int>
# 1 Coleoptera    103
# 2 Diptera        76
# 3 Hemiptera      68
# 4 Hymenoptera    54
# 5 Lepidoptera    60


pred_data <- predicts

##### For each family, remove it from the dataset and rerun the analysis ####

# load in the predicts dataset including the climate data so it can be merged
pred_ori <- readRDS("4_PREDICTSMatchClimateIndex/PREDICTSSites_Climate.rds")
pred_ori <- pred_ori@data


#### start for loop ####

plot_data_ab <- NULL
plot_data_sr <- NULL

for(family in fams){
  
predicts <- pred_data[!pred_data$Family == fams[family], ]  


#### 2. Calculate site metrics and prepare dataset ####

# Split predicts into separate data frames according to insect Order 

# use split function to split the predicts data frame into 6 data frames (1/Order)
OrderName <- paste0("",predicts$Order)

by_Order <- split(predicts,OrderName)

# extract data frames from list into global environment
list2env(by_Order,globalenv())

# Calculate site metrics of diversity for each order, include extra columns:
Coleoptera <- droplevels(Coleoptera)
Coleoptera <- SiteMetrics(diversity = Coleoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion"))
Diptera <- droplevels(Diptera)
Diptera <- SiteMetrics(diversity = Diptera,
                       extra.cols = c("Predominant_land_use",
                                      "SSB","SSBS", "Biome","Order","UN_subregion"))
Hemiptera <- droplevels(Hemiptera)
Hemiptera <- SiteMetrics(diversity = Hemiptera,
                         extra.cols = c("Predominant_land_use",
                                        "SSB","SSBS", "Biome","Order","UN_subregion"))
Hymenoptera <- droplevels(Hymenoptera)
Hymenoptera <- SiteMetrics(diversity = Hymenoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order","UN_subregion"))
Lepidoptera <- droplevels(Lepidoptera)
Lepidoptera <- SiteMetrics(diversity = Lepidoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order","UN_subregion"))

# merge all sites_Order data frames into one called "sites"
# merge using rbind()
sites <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera)

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


# remove the urban sites and sites that are NA in LUI
sites <- sites[!sites$LUI == "Urban", ]
sites <- sites[!is.na(sites$LUI), ]


sites <- droplevels(sites)

# transform abundance values 
sites$LogAbund <- log(sites$Total_abundance+1)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]


# create a new variable designating sites as Tropical or Non-tropical
# assign new variable for tropical/temperate, convert to factor, and filter out NA
sites$Realm <- ifelse(sites$Latitude >= -23.5 & sites$Latitude <= 23.5, "Tropical", "NonTropical")
sites$Realm <- factor(sites$Realm, levels = c("NonTropical", "Tropical"))
sites <- sites %>%
  filter(!is.na(Realm))

#### 3. Merge sites with main dataset to add the climate variables ####

pred_ori2 <- pred_ori[, c(15, 29:35)]
pred_ori2 <- unique(pred_ori2)

sites <- left_join(sites, pred_ori2, by = "SSBS")


# set LUI as factor and set reference level
sites$LUI <- factor(sites$LUI)
sites$LUI <- relevel(sites$LUI,ref="Primary vegetation")

# rescale the climate variables
sites$StdTmeanAnomalyRS <- StdCenterPredictor(sites$StdTmeanAnomaly)
sites$StdTmaxAnomalyRS <- StdCenterPredictor(sites$StdTmaxAnomaly)

# rescaling abundance and log values
# CO note changed from RescaleAbundance to RescaleAbundance2 to scale within Study AND Order
sites <- RescaleAbundance2(sites)

# Charlie added this line as later bits were throwing errors
sites <- droplevels(sites)

# some of the climate values are NA since they do not meet the thresholds
sites <- sites[!is.na(sites$avg_temp), ]



#### 4. rerun the analysis with subsetted data ####

# i. Abundance, mean anomaly

model_data <- sites[!is.na(sites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))


# ii. Richness, mean anomaly

model_data <- sites[!is.na(sites$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS"))


#### 5. plotting data ####

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

##### a. Abundance, Mean Anomaly ####

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = sites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))


# set quantiles by Order
# coleoptera
C_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Coleoptera"],
  probs = exclQuantiles)

# diptera
D_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Diptera"],
  probs = exclQuantiles)

# hemiptera
He_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Hemiptera"],
  probs = exclQuantiles)

# hymenoptera
Hy_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Hymenoptera"],
  probs = exclQuantiles)

# lepidoptera
L_QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High" & MeanAnomalyModelAbund$data$Order == "Lepidoptera"],
  probs = exclQuantiles)


# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# create list of matrices
number_of_chunks = 5
list_a.preds.tmean <- lapply(seq(1, NROW(a.preds.tmean), ceiling(NROW(a.preds.tmean)/number_of_chunks)),
                             function(i) a.preds.tmean[i:min(i + ceiling(NROW(a.preds.tmean)/number_of_chunks) - 1, NROW(a.preds.tmean)),])

# name them
names(list_a.preds.tmean) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")

# sweep out according to the refRow
list_a.preds.tmean <- lapply(list_a.preds.tmean,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[refRow[1],],FUN="/") 
})


# print to global environment
list2env(list_a.preds.tmean,globalenv())

# split nd by order
Order<- paste0("nd_",nd$Order)
# create a list of data frames
by_Order <- split(nd,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QPV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QPV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QSV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QSV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QAL[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QAL[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS < C_QAH[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS > C_QAH[2])),] <- NA

Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < D_QPV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > D_QPV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < D_QSV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > D_QSV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS < D_QAL[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS > D_QAL[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS < D_QAH[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS > D_QAH[2])),] <- NA

Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QPV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QPV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QSV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QSV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QAL[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QAL[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS < He_QAH[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS > He_QAH[2])),] <- NA

Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QPV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QPV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QSV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QSV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QAL[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QAL[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS < Hy_QAH[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS > Hy_QAH[2])),] <- NA

Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS < L_QPV[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS > L_QPV[2])),] <- NA
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
nd_Coleoptera$LUI <- factor(nd_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Diptera$LUI <- factor(nd_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hemiptera$LUI <- factor(nd_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hymenoptera$LUI <- factor(nd_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Lepidoptera$LUI <- factor(nd_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


all_groups <- rbind(nd_Coleoptera, nd_Diptera, nd_Hemiptera, nd_Hymenoptera, nd_Lepidoptera)

all_groups$family <- family

plot_data_ab <- rbind(plot_data_ab, all_groups)



#### 4. Richness, Mean Anomaly ####

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = sites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$LUI=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# set quantiles
# coleoptera
C_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)
C_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Coleoptera"],
  probs = exclQuantiles)

# diptera
D_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)
D_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Diptera"],
  probs = exclQuantiles)

# hemiptera
He_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)
He_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Hemiptera"],
  probs = exclQuantiles)

# hymenoptera
Hy_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)
Hy_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Hymenoptera"],
  probs = exclQuantiles)

# lepidoptera
L_QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)
L_QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High" & MeanAnomalyModelRich$data$Order == "Lepidoptera"],
  probs = exclQuantiles)


# predict the results
sr.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd2)

# back transform the abundance values
sr.preds.tmean <- exp(sr.preds.tmean)-0.01

# create list of matrices
number_of_chunks = 5
list_sr.preds.tmean <- lapply(seq(1, NROW(sr.preds.tmean), ceiling(NROW(sr.preds.tmean)/number_of_chunks)),
                              function(i) sr.preds.tmean[i:min(i + ceiling(NROW(sr.preds.tmean)/number_of_chunks) - 1, NROW(sr.preds.tmean)),])
# name them
names(list_sr.preds.tmean) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")

# sweep out according to the refRow
list_sr.preds.tmean <- lapply(list_sr.preds.tmean,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[refRow[1],],FUN="/") 
})

list2env(list_sr.preds.tmean,globalenv())

# split nd by order
Order<- paste0("nd2_",nd2$Order)
# create a list of data frames
by_Order <- split(nd2,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QPV[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QPV[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QSV[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QSV[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QAL[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QAL[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS < C_QAH[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS > C_QAH[2])),] <- NA

Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < D_QPV[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > D_QPV[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < D_QSV[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > D_QSV[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS < D_QAL[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS > D_QAL[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS < D_QAH[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS > D_QAH[2])),] <- NA

Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QPV[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QPV[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QSV[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QSV[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QAL[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QAL[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS < He_QAH[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS > He_QAH[2])),] <- NA

Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QPV[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QPV[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QSV[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QSV[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QAL[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QAL[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS < Hy_QAH[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS > Hy_QAH[2])),] <- NA

Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QPV[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QPV[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QSV[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QSV[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QAL[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QAL[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS < L_QAH[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS > L_QAH[2])),] <- NA



# Get the median, upper and lower quants for the plot
nd2_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd2_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                  FUN = median,na.rm=TRUE))*100)-100
nd2_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd2_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd2_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd2_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd2_Coleoptera$LUI <- factor(nd2_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Diptera$LUI <- factor(nd2_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Hemiptera$LUI <- factor(nd2_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Hymenoptera$LUI <- factor(nd2_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Lepidoptera$LUI <- factor(nd2_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


all_groups <- rbind(nd2_Coleoptera, nd2_Diptera, nd2_Hemiptera, nd2_Hymenoptera, nd2_Lepidoptera)

all_groups$family <- family

plot_data_sr <- rbind(plot_data_sr, all_groups)


}


# plot

p_coleoptera <- ggplot(data = nd_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Coleoptera$PredLower, ymax = nd_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") + 
  ggtitle("Coleoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 


p_diptera <- ggplot(data = nd_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Diptera$PredLower, ymax = nd_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Diptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 


p_hemiptera <- ggplot(data = nd_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hemiptera$PredLower, ymax = nd_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hemiptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))


p_hymenoptera <- ggplot(data = nd_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hymenoptera$PredLower, ymax = nd_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hymenoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_lepidoptera <- ggplot(data = nd_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Lepidoptera$PredLower, ymax = nd_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Lepidoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 


# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 7), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomAbund <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,legend, ncol=3)

# save plot (pdf)
ggsave(filename = paste0(plotDir, "MeanAnomAbund.pdf"), plot = MeanAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("MeanAnomAbund.jpeg", device ="jpeg", path = plotDir, width=20, height=15, units="cm", dpi = 350)





# plot
p_coleoptera <- ggplot(data = nd2_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Coleoptera$PredLower, ymax = nd2_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Coleoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 

p_diptera <- ggplot(data = nd2_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Diptera$PredLower, ymax = nd2_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Diptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_hemiptera <- ggplot(data = nd2_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hemiptera$PredLower, ymax = nd2_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hemiptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_hymenoptera <- ggplot(data = nd2_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hymenoptera$PredLower, ymax = nd2_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Hymenoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25))

p_lepidoptera <- ggplot(data = nd2_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Lepidoptera$PredLower, ymax = nd2_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.25) +
  scale_fill_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual('Land-use', values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  ggtitle("Lepidoptera") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 8, face = 'bold'),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25)) 

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 7), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomRich <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,legend, ncol=3)

# save plot (pdf)
ggsave(filename = paste0(plotDir, "MeanAnomRich.pdf"), plot = MeanAnomRich, width = 200, height = 150, units = "mm", dpi = 300)

# save plot (jpeg)
ggsave("MeanAnomRich.jpeg", device ="jpeg", path = plotDir, width=20, height=15, units="cm", dpi = 350)

