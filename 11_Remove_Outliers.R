##%######################################################%##
#                                                          #
####        11. Additional test: Remove outliers        ####
#                                                          #
##%######################################################%##


# This script checks for outliers in the models, removes these from the 
# dataset and reruns the models to see if they have an effect on the results.


# Here looking for influential studies using the influence.ME package

rm(list = ls())

# directories
moddir <- "5_RunLUIClimateModels/"
outdir <- "10_Additional_Tests/Outliers/"
if(!dir.exists(outdir)) dir.create(outdir)

sink(paste0(outdir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# libraries
library(performance)
source("0_Functions.R")
library(influence.ME)
library(StatisticalModels)


# load in models
load(paste0(moddir, "MeanAnomalyModelAbund.rdata")) # MeanAnomalyModelAbund
load(paste0(moddir, "MeanAnomalyModelRich.rdata"))  # MeanAnomalyModelRich


############################################################
#                                                          #
####                 Assess outliers                    ####
#                                                          #
############################################################


# 1. abundance mean anomaly

# result1 <- cooks.distance(MeanAnomalyModelAbund$model, sort = T)
# range(result1)
# summary(result1 < 1)# all true
# result1[result1 > 0.4]

# extract the data used to run the model
modelData <- MeanAnomalyModelAbund$data

# number of unique studies
length(unique(modelData$SS)) # 235

# rerun the model using lme4 to avoid errors when using the influence function
mod1 <- lmer(formula = LogAbund ~ LUI * StdTmeanAnomalyRS * Order + (1 | SS) + (1 | SSB), modelData)



alt.est.a <- influence(mod1, "SS")
plot(alt.est.a, which = "cook", sort = T) # takes a long time to run
result1 <- cooks.distance(alt.est.a, sort = T)

# studies with cook's distance >= 1
# SC1_2005__Richardson 
# HW1_2011__Summerville 
# AD1_2008__Billeter 



# 2. Richness mean anomaly

modelData <- MeanAnomalyModelRich$data

length(unique(modelData$SS)) # 263

# rerun the model using lme4 to avoid errors when using the influence function
mod2 <- glmer(formula = Species_richness~LUI * StdTmeanAnomalyRS * Order+(1|SS)+(1|SSB)+(1|SSBS), modelData, family = "poisson")


# additional info needed for the influence function
# maxIters=10000
# optimizer="bobyqa"
# fitFamily = "poisson"
alt.est.a <- influence(mod2, "SS") # takes a long time to run
plot(alt.est.a, which = "cook", sort = T)
result2 <- cooks.distance(alt.est.a, sort = T)






##%######################################################%##
#                                                          #
####      Rerun the models removing these outliers      ####
#                                                          #
##%######################################################%##

datadir <- "6_RunLUClimateModels/"

predictsSites <- readRDS(file = paste0(datadir,"PREDICTSSiteData.rds"))



# 1. Abundance, mean anomaly

# list of studies that had larger Cook's distance values
outliers <- c("SC1_2011__Meijer 1", "CC1_2007__Ewers 1")

# subset the data to exclude the potential outliers
model_data <- predictsSites[!predictsSites$SS %in% outliers, ] # 5710 rows

model_data <- model_data[!is.na(model_data$LogAbund), ] #5374 rows

length(unique(model_data$SS)) # 242

MeanAnomalyModelAbund <- GLMERSelect(modelData = model_data,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outdir, "/MeanAnomalyModelAbund_rmout.rdata"))
#load(file = paste0(outdir, "/MeanAnomalyModelAbund_rmout.rdata"))



# 2. Richness, mean anomaly


# list of studies that had larger Cook's distance values
outliers <- c("SC1_2011__Meijer 1", "CC1_2007__Ewers 1", "HP1_2014__Gray 1")


# subset the data to exclude the potential outliers
model_data <- predictsSites[!predictsSites$SS %in% outliers, ] # 5419 rows

model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]

nrow(predictsSites[is.na(predictsSites$StdTmeanAnomalyRS), ])

length(unique(model_data$SS)) # 260

MeanAnomalyModelRich <- GLMERSelect(modelData = model_data,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("SSBS"))

# save model output
save(MeanAnomalyModelRich, file = paste0(outdir, "/MeanAnomalyModelRich_rmout.rdata"))
#load(file = paste0(outdir, "/MeanAnomalyModelRich_rmout.rdata"))

MeanAnomalyModelRich$model


##%######################################################%##
#                                                          #
####                   Replot Figures                   ####
#                                                          #
##%######################################################%##


## copy over code from other script





t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()

