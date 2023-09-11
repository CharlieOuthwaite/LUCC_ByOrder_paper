############################################################
#                                                          #
#            Models without rescaling abundance            #
#                                                          #
############################################################

# In this script, the models for land use and land use/ climate impacts
# are rerun for total abundance, however abundance values are not rescaled. 


# ensure working directory is clear
rm(list = ls())

# set up directories
inDir <- "1_CheckPrepareData/"
outDir <- "10_Additional_Tests/Without_rescaling_abundance/"
predsDir <- "7_Predictions/"
dataDir <- "5_RunLUIClimateModels/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(predsDir)) dir.create(predsDir)

# load libraries
packages_model <- c("StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS","webshot")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")



############################################################
#                                                          #
####               LUI and STA models                   ####
#                                                          #
############################################################


# load in the data
predictsSites <- readRDS(file = paste0(dataDir,"PREDICTSSitesClimate_Data.rds"))

# create new variable where abundance data is not rescaled
predictsSites$LogAbund_norescaling <- log(predictsSites$Total_abundance + 0.01)

# hist(predictsSites$LogAbund_norescaling)
# hist(predictsSites$LogAbund)


# i. Abundance, mean anomaly including interaction

model_data <- predictsSites[!is.na(predictsSites$LogAbund_norescaling), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund_noresc <- GLMER(modelData = model_data,responseVar = "LogAbund_norescaling",fitFamily = "gaussian",
                               fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund_noresc$model)


# save the model output
save(MeanAnomalyModelAbund_noresc, file = paste0(outDir, "MeanAnomalyModelAbund_incOrder_norescaling.rdata"))



# i. Abundance, mean anomaly excluding interaction with order

MeanAnomalyModelAbund2_noresc <- GLMER(modelData = model_data,responseVar = "LogAbund_norescaling",fitFamily = "gaussian",
                                fixedStruct = "LUI * StdTmeanAnomalyRS",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund2_noresc$model)


# save the model output
save(MeanAnomalyModelAbund2_noresc, file = paste0(outDir, "MeanAnomalyModelAbund_noOrder_norescaling.rdata"))



############################################################
#                                                          #
####              Abundance predictions                 ####   
#                                                          #
############################################################


# what is the rescaled value of STA of 1
BackTransformCentreredPredictor(transformedX = 0.999, originalX = predictsSites$StdTmeanAnomaly) # 0.999 gives about 1 

# what is the rescaled value of STA of 0
BackTransformCentreredPredictor(transformedX = -1.39, originalX = predictsSites$StdTmeanAnomaly) # -1.39 gives about 0 


#### 1. No rescaling but including interaction with Order ####

# set up table for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.999),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund_noresc$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund_norescaling <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
# row for each order

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_noresc$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to dataframe
result.ab <- as.data.frame(a.preds.tmean)

# add in the LU info
result.ab$LUI <- nd$LUI

# add in order info
result.ab$Order <- nd$Order

# express as a percentage of primary
# break into Orders
Order<- paste0("",result.ab$Order)
list.result.ab <- split(result.ab,Order)
list2env(list.result.ab,globalenv())

# convert to percentage difference from primary vegetation
Coleoptera <- as.matrix(Coleoptera[, 1:1000])
col.preds <- sweep(x = Coleoptera, MARGIN = 2, STATS = Coleoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = col.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

colres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
colres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
colres$Order <- "Coleoptera"

# convert to percentage difference from primary vegetation
Diptera <- as.matrix(Diptera[, 1:1000])
dip.preds <- sweep(x = Diptera, MARGIN = 2, STATS = Diptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = dip.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

dipres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
dipres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
dipres$Order <- "Diptera"

# convert to percentage difference from primary vegetation
Hemiptera <- as.matrix(Hemiptera[, 1:1000])
hem.preds <- sweep(x = Hemiptera, MARGIN = 2, STATS = Hemiptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hem.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hemres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hemres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
hemres$Order <- "Hemiptera"

# convert to percentage difference from primary vegetation
Hymenoptera <- as.matrix(Hymenoptera[, 1:1000])
hym.preds <- sweep(x = Hymenoptera, MARGIN = 2, STATS = Hymenoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hym.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hymres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hymres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
hymres$Order <- "Hymenoptera"

# convert to percentage difference from primary vegetation
Lepidoptera <- as.matrix(Lepidoptera[, 1:1000])
lep.preds <- sweep(x = Lepidoptera, MARGIN = 2, STATS = Lepidoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = lep.preds,MARGIN = 1, FUN = median))*100)-100
grp.upper <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.975))*100)-100
grp.lower <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.025))*100)-100

lepres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
lepres$LUI <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_Low", "Agriculture_High", "Agriculture_High")
lepres$Order <- "Lepidoptera"

# put it back together
result.ab <- rbind(colres, dipres, hemres, hymres, lepres)
result.ab$Metric <- "Total abundance"

all_res <- result.ab

#### 2. No rescaling but excluding interaction with Order ####


# reference is primary with 0 climate change so have 0 for that row
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.999),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund2_noresc$data$LUI)))

# set values for prediction
nd$LogAbund_norescaling <- 0

# back transform the predictors and round
nd$StdTmeanAnomaly <- round(BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly))

# predict the results
result.ab2 <- PredictGLMERRandIter(model = MeanAnomalyModelAbund2_noresc$model, data = nd)


# backtransform
result.ab2 <- exp(result.ab2)-0.01

# convert to percentage difference from primary vegetation
result.ab2 <- sweep(x = result.ab2, MARGIN = 2, STATS = result.ab2[1,], FUN = '/')

# get quantiles
a.preds.median <- ((apply(X = result.ab2,MARGIN = 1,FUN = median))*100)-100
a.preds.upper <- ((apply(X = result.ab2,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
a.preds.lower <- ((apply(X = result.ab2,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


# combine data into one table for plotting
abun_res <- as.data.frame(cbind(a.preds.median, a.preds.lower, a.preds.upper))
colnames(abun_res) <- c("median", "lower", "upper")
abun_res$Metric <- "Total abundance"
abun_res$LUI <- factor(c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low","Agriculture_Low", "Agriculture_High",  "Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))

abun_res$STA <- nd$StdTmeanAnomaly

abun_res$Study <- "This study"
abun_res$Order <- "All insects"
abun_res$Fixed_effs <- "Land use and climate"
abun_res$Realm <- "Global"

res <- abun_res[, c(7, 4, 8:10, 5, 1:3, 6)]

all_res$Study <- "This study"
all_res$Fixed_effs <- "Land use and climate"
all_res$Realm <- "Global"
all_res$STA <- c(0, 1)
names(all_res)[1:3] <- c("median", "lower", "upper")
all_res <- all_res[, c(7, 6, 5, 8, 9, 4, 1:3, 10)]

all_res <- rbind(all_res, res)

# save table
write.csv(all_res, file = paste0(outDir, "/percentage_change_LU_CC_abundance_norescaling.csv"), row.names = F)



all_res$Order <- factor(all_res$Order, levels = c("All insects", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"))

all_res$STA <- factor(all_res$STA, levels = c("0", "1"))

all_res$LUI <- factor(all_res$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


# create point and error bar plot
ggplot(data = all_res, aes(col = LUI, group = STA)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  geom_point(aes(x = LUI, y = median, shape = STA), size = 1.5, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = LUI, ymin = lower, ymax = upper), position= position_dodge(width = 1), size = 0.5, width = 0.2)+
  facet_wrap(~ Order) +
  xlab("") +
  scale_y_continuous(limits = c(-100, 180), breaks = scales::pretty_breaks(n = 10)) +
  ylab("Percentage change in total abundance (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00"), guide = "none") +
  scale_shape_manual(values=c(16, 17, 18, 15, 0, 1), name = "STA")+
  theme(legend.position = "bottom", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
        axis.title = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2), 
        text = element_text(size = 8), 
        legend.key=element_blank(), 
        #legend.title = element_blank(),
        legend.text = element_text(size = 8))


ggsave(filename = paste0(outDir, "Comparison_LUSTA_Abun_norescaling.pdf"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "Comparison_LUSTA_Abun_norescaling.jpeg"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 600)



