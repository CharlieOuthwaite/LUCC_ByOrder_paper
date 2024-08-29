############################################################
#                                                          #
#           Predictions from both sets of models           #
#                                                          #
############################################################


# In this script, I predict the biodiversity change for different combinations of
# land use and STA for both the models from Outhwaite et al (2022) and those
# tested here (including and excluding the order-level interactions)

# clear environment
rm(list = ls())

# directories
predictsDataDir <- "5_RunLUIClimateModels/"
moddir1 <- "2_RunSimpleLUIModel/"
moddir2 <- "5_RunLUIClimateModels/"
outDir <- "7_Predictions/Model_Comparisons/"
if(!dir.exists(outDir)) dir.create(outDir)

# load libraries
library(StatisticalModels)
library(predictsFunctions)
library(ggplot2)
library(paletteer) 
source('0_Functions.R')

# read in the predicts data
predictsSites <- readRDS(paste0(predictsDataDir, "PREDICTSSitesClimate_Data.rds")) # 8858 rows


##%######################################################%##
#                                                          #
####            Hyp 1: land use effect only             ####
#                                                          #
##%######################################################%##

#### Models from this study excluding the interaction with order ####

load(file = paste0(moddir1, "Abundance_landuse_model_noOrder.rdata")) # am3
load(file = paste0(moddir1, "Richness_landuse_model_noOrder.rdata")) # sm3


# create dataframe for values to predict response to
nd <- data.frame(LUI=factor(c("Primary vegetation","Secondary vegetation",
                              "Agriculture_Low","Agriculture_High"),
                            levels = levels(sm3$data$LUI)),
                 Species_richness=0,
                 LogAbund=0)

## species richness predictions ##
s.preds <- PredictGLMERRandIter(model = sm3$model, data = nd)

s.preds <- exp(s.preds)

# convert to percentage difference from primary vegetation
s.preds <- sweep(x = s.preds, MARGIN = 2, STATS = s.preds[1,], FUN = '/')

# get quantiles
s.preds.median <- ((apply(X = s.preds,MARGIN = 1,FUN = median))*100)-100
s.preds.upper <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
s.preds.lower <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


## abundance predictions ##

a.preds <- PredictGLMERRandIter(model = am3$model,data = nd)

a.preds <- exp(a.preds)-0.01

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
abun_res$metric <- "total abundance"
rich_res$metric <- "species richness"
abun_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))
rich_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))

abun_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA
rich_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA

# combine results
noOrder_res <- rbind(abun_res, rich_res)

##%######################################################%##
#                                                          #
####  plot species richness and abundance predictions   ####
#                                                          #
##%######################################################%##


# point plots

p1 <- ggplot(data = abun_res) +
  geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
  geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  xlab("") +
  ylab("Change in total abundance (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
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
        axis.ticks = element_line(linewidth = 0.2), 
        axis.line = element_line(linewidth = 0.2))


p2 <- ggplot(data = rich_res) +
  geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
  geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  xlab("") +
  ylab("Change in species richness (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
  ylim(c(-50, 0)) +
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
        axis.ticks = element_line(linewidth = 0.2), 
        axis.line = element_line(linewidth = 0.2))


p3 <- cowplot::plot_grid(p1, p2)

ggsave(filename = paste0(outDir, "FIGURE_landuse_noOrder.pdf"), plot = last_plot(), width = 120, height = 60, units = "mm", dpi = 300)


##### Models from this study including interaction with order ####

# load model outputs
load(file = paste0(moddir1, "Abundance_landuse_model.rdata")) # am3.3
load(file = paste0(moddir1, "Richness_landuse_model.rdata")) # sm3.3

# create table for predictions
data_tab <- expand.grid(LUI = factor(c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"), levels = levels(am3.3$data$LUI)), 
                       Order = factor(c("Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"), levels = levels(am3.3$data$Order)),
                       LogAbund = 0,
                       Species_richness = 0)


# predict the results
result.sr <- PredictGLMERRandIter(model = sm3.3$model, data = data_tab)

# backtransform
result.sr <- exp(result.sr)

# convert to dataframe
result.sr <- as.data.frame(result.sr)

# add in the LUI info
result.sr$LUI <- data_tab$LUI

# add in the Order info
result.sr$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr$Order)
list.result.sr <- split(result.sr,Order)
list2env(list.result.sr,globalenv())

# convert to percentage difference from primary vegetation
Coleoptera <- as.matrix(Coleoptera[, 1:1000])
col.preds <- sweep(x = Coleoptera, MARGIN = 2, STATS = Coleoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = col.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = col.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

colres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
colres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
colres$Order <- "Coleoptera"

# convert to percentage difference from primary vegetation
Diptera <- as.matrix(Diptera[, 1:1000])
dip.preds <- sweep(x = Diptera, MARGIN = 2, STATS = Diptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = dip.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

dipres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
dipres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
dipres$Order <- "Diptera"

# convert to percentage difference from primary vegetation
Hemiptera <- as.matrix(Hemiptera[, 1:1000])
hem.preds <- sweep(x = Hemiptera, MARGIN = 2, STATS = Hemiptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hem.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hemres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hemres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hemres$Order <- "Hemiptera"

# convert to percentage difference from primary vegetation
Hymenoptera <- as.matrix(Hymenoptera[, 1:1000])
hym.preds <- sweep(x = Hymenoptera, MARGIN = 2, STATS = Hymenoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hym.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hymres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hymres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hymres$Order <- "Hymenoptera"

# convert to percentage difference from primary vegetation
Lepidoptera <- as.matrix(Lepidoptera[, 1:1000])
lep.preds <- sweep(x = Lepidoptera, MARGIN = 2, STATS = Lepidoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = lep.preds,MARGIN = 1, FUN = median))*100)-100
grp.upper <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.975))*100)-100
grp.lower <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.025))*100)-100

lepres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
lepres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
lepres$Order <- "Lepidoptera"

# put it back together
result.sr <- rbind(colres, dipres, hemres, hymres, lepres)
result.sr$metric <- "species richness"

# now for the abundance model  
result.ab <- PredictGLMERRandIter(model = am3.3$model, data = data_tab)

# backtransform
result.ab <- exp(result.ab)-0.01

# convert to dataframe
result.ab <- as.data.frame(result.ab)

# add in the LU info
result.ab$LUI <- data_tab$LUI

# add in order info
result.ab$Order <- data_tab$Order

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
colres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
colres$Order <- "Coleoptera"

# convert to percentage difference from primary vegetation
Diptera <- as.matrix(Diptera[, 1:1000])
dip.preds <- sweep(x = Diptera, MARGIN = 2, STATS = Diptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = dip.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = dip.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

dipres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
dipres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
dipres$Order <- "Diptera"

# convert to percentage difference from primary vegetation
Hemiptera <- as.matrix(Hemiptera[, 1:1000])
hem.preds <- sweep(x = Hemiptera, MARGIN = 2, STATS = Hemiptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hem.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hem.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hemres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hemres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hemres$Order <- "Hemiptera"

# convert to percentage difference from primary vegetation
Hymenoptera <- as.matrix(Hymenoptera[, 1:1000])
hym.preds <- sweep(x = Hymenoptera, MARGIN = 2, STATS = Hymenoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = hym.preds,MARGIN = 1,FUN = median))*100)-100
grp.upper <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
grp.lower <- ((apply(X = hym.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

hymres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
hymres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
hymres$Order <- "Hymenoptera"

# convert to percentage difference from primary vegetation
Lepidoptera <- as.matrix(Lepidoptera[, 1:1000])
lep.preds <- sweep(x = Lepidoptera, MARGIN = 2, STATS = Lepidoptera[1,], FUN = '/')

# get quantiles
grp.median <- ((apply(X = lep.preds,MARGIN = 1, FUN = median))*100)-100
grp.upper <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.975))*100)-100
grp.lower <- ((apply(X = lep.preds,MARGIN = 1, FUN = quantile, probs = 0.025))*100)-100

lepres <- as.data.frame(cbind(grp.median, grp.upper, grp.lower))
lepres$LUI <- c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High")
lepres$Order <- "Lepidoptera"

# put it back together
result.ab <- rbind(colres, dipres, hemres, hymres, lepres)
result.ab$metric <- "total abundance"

all_res <- rbind(result.ab, result.sr)


# match up noOrder results for combining
noOrder_res$Order <- "All insects"
noOrder_res <- noOrder_res[, c(1:3, 5, 6, 4)]
names(noOrder_res) <- c("grp.median", "grp.lower", "grp.upper", "Metric", "LUI", "Order") 

## add in results from noOrder models
all_res <- rbind(all_res, noOrder_res)

# replace names
all_res$LUI <- sub("Agriculture_Low", "Low-intensity agriculture", all_res$LUI)
all_res$LUI <- sub("Agriculture_High", "High-intensity agriculture", all_res$LUI)

all_res$Study <- "This study"
all_res$Fixed_effs <- "Land use only"

# organise columns

all_res <- all_res[, c(7, 6, 5, 8, 4, 1, 3, 2)]

# rename columns
names(all_res) <- c("Study", "Metric", "Order", "Fixed_effs", "LUI", "Median", "Lower_CI", "Upper_CI")

# round percentages to 3 dp
all_res[, 6:8] <- round(all_res[, 6:8], 3)

# save table
write.csv(all_res, file = paste0(outDir, "/TABLE_percentage_change_LU_Order.csv"), row.names = F)


##%######################################################%##
#                                                          #
####        Figure for percentage changes by LU         ####
#                                                          #
##%######################################################%##

all_res$LUI <- factor(all_res$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Low-intensity agriculture", "High-intensity agriculture"))
all_res$Order <- factor(all_res$Order, levels = c("All insects", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"))

all_res$Metric <- sub("species richness", "Species richness", all_res$Metric)
all_res$Metric <- sub("total abundance", "Total abundance", all_res$Metric)

# create point and error bar plot
ggplot(data = all_res, aes(col = LUI, group = LUI)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  geom_point(aes(x = Order, y = Median, col = LUI), size = 2, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = Order, ymin = Lower_CI, ymax = Upper_CI), position= position_dodge(width = 1), size = 0.5, width = 0.2)+
  facet_wrap(~ Metric) +
  xlab("") +
  scale_y_continuous(limits = c(-80, 80), breaks = scales::pretty_breaks(n = 10)) +
  ylab("Percentage change (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
  theme(legend.position = "bottom", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2), 
        text = element_text(size = 10), 
        legend.key=element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10))


ggsave(filename = paste0(outDir, "FIGURE_2_Comparison_LU_only_incNoOrder.pdf"), plot = last_plot(), width = 200, height = 130, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "FIGURE_2_Comparison_LU_only_incNoOrder.jpeg"), plot = last_plot(), width = 200, height = 130, units = "mm", dpi = 300)



##%######################################################%##
#                                                          #
####  Hyp 2: Land use and climate anomaly interaction   ####
#                                                          #
##%######################################################%##


#### This study including interaction with Order ####
# need predictions for each land use for 0 and 1 STA for each order

# load in models
load(file = paste0(moddir2, "MeanAnomalyModelAbund.rdata")) # MeanAnomalyModelAbund
load(file = paste0(moddir2, "MeanAnomalyModelRich.rdata")) # MeanAnomlayModelRich

#### abundance predictions ####

# what is the rescaled value of STA of 1
BackTransformCentreredPredictor(transformedX = 0.999, originalX = predictsSites$StdTmeanAnomaly) # 0.999 gives about 1 

# what is the rescaled value of STA of 0
BackTransformCentreredPredictor(transformedX = -1.39, originalX = predictsSites$StdTmeanAnomaly) # -1.39 gives about 0 

# set up table for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.999),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
# row for each order

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)

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

#### species richness predictions ####

# set up table for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.999),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
# row for each order, its the 1st row each time

# predict the results
s.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd)

# back transform the abundance values
s.preds.tmean <- exp(s.preds.tmean)

# convert to dataframe
result.sr <- as.data.frame(s.preds.tmean)

# add in the LU info
result.sr$LUI <- nd$LUI

# add in order info
result.sr$Order <- nd$Order

# express as a percentage of primary
# break into Orders
Order<- paste0("",result.sr$Order)
list.result.sr <- split(result.sr,Order)
list2env(list.result.sr,globalenv())

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
result.sr <- rbind(colres, dipres, hemres, hymres, lepres)
result.sr$Metric <- "Species Richness"


all_res <- rbind(result.ab, result.sr)

all_res$Study <- "This study"

all_res$Fixed_effs <- "Land use and climate"

all_res$Realm <- "Global"

all_res$STA <- c(0, 1)

all_res <- all_res[ , c(7,6,5, 8, 9, 4, 1, 3, 2, 10)]

names(all_res)[7:9] <- c("Median", "Lower_CI", "Upper_CI")

#### This study excluding interaction with Order ####

# load in models
load(file = paste0(moddir2, "/MeanAnomalyModelAbund_noOrder.rdata")) # MeanAnomalyModelAbund2
load(file = paste0(moddir2, "/MeanAnomalyModelRich_noOrder.rdata")) # MeanAnomalyModelRich2


# create matrix for predictions
# Primary, Low, High
# SCA = 0, 1
# abun and richness = 0

# what is the rescaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = 0.999, originalX = predictsSites$StdTmeanAnomaly) # 0.999 gives about 1 

# what is the rescaled value of SCA of 0
BackTransformCentreredPredictor(transformedX = -1.39, originalX = predictsSites$StdTmeanAnomaly) # -1.39 gives about 0 

# reference is primary with 0 climate change so have 0 for that row
nd <- expand.grid(
  StdTmeanAnomalyRS= c(-1.39, 0.999),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich2$data$LUI)))

# set values for prediction
nd$Species_richness <- 0
nd$LogAbund <- 0

# back transform the predictors and round
nd$StdTmeanAnomaly <- round(BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly))

# predict the results
result.sr <- PredictGLMERRandIter(model = MeanAnomalyModelRich2$model, data = nd)

# backtransform
result.sr <- exp(result.sr)

# convert to percentage difference from primary vegetation
result.sr <- sweep(x = result.sr, MARGIN = 2, STATS = result.sr[1,], FUN = '/')

# get quantiles
s.preds.median <- ((apply(X = result.sr,MARGIN = 1,FUN = median))*100)-100
s.preds.upper <- ((apply(X = result.sr,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
s.preds.lower <- ((apply(X = result.sr,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


# now for the abundance model  
result.ab <- PredictGLMERRandIter(model = MeanAnomalyModelAbund2$model, data = nd)

# backtransform
result.ab <- exp(result.ab)-0.01

# convert to percentage difference from primary vegetation
result.ab <- sweep(x = result.ab, MARGIN = 2, STATS = result.ab[1,], FUN = '/')

# get quantiles
a.preds.median <- ((apply(X = result.ab,MARGIN = 1,FUN = median))*100)-100
a.preds.upper <- ((apply(X = result.ab,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
a.preds.lower <- ((apply(X = result.ab,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


# combine data into one table for plotting
abun_res <- as.data.frame(cbind(a.preds.median, a.preds.lower, a.preds.upper))
rich_res <- as.data.frame(cbind(s.preds.median, s.preds.lower, s.preds.upper))
colnames(abun_res) <- c("median", "lower", "upper")
colnames(rich_res) <- c("median", "lower", "upper")
abun_res$metric <- "Total abundance"
rich_res$metric <- "Species richness"
abun_res$LU <- factor(c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low","Agriculture_Low", "Agriculture_High",  "Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))
rich_res$LU <- factor(c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation", "Agriculture_Low","Agriculture_Low", "Agriculture_High",  "Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))

abun_res$STA <- nd$StdTmeanAnomaly
rich_res$STA <- nd$StdTmeanAnomaly

# add info and match columns to all_res
res <- rbind(abun_res, rich_res)

res$Study <- "This study"
res$Order <- "All insects"
res$Fixed_effs <- "Land use and climate"
res$Realm <- "Global"

res <- res[, c(7, 4, 8:10, 5, 1:3, 6)]
names(res) <- names(all_res)

all_res <- rbind(all_res, res)

# save table
write.csv(all_res, file = paste0(outDir, "/TABLE_percentage_change_LU_CC_incNoOrder.csv"), row.names = F)


##%######################################################%##
#                                                          #
####     Figure for percentage changes by LU and STA    ####
#                                                          #
##%######################################################%##

final_res <- all_res

final_res$Order <- factor(final_res$Order, levels = c("All insects", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"))

final_res$STA <- factor(final_res$STA, levels = c("0", "1"))

final_res$LUI <- factor(final_res$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

plot_data <- final_res[final_res$Metric == "Total abundance", ]


# create point and error bar plot
ggplot(data = plot_data, aes(col = LUI, group = STA)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  geom_point(aes(x = LUI, y = Median, shape = STA), size = 1.5, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = LUI, ymin = Lower_CI, ymax = Upper_CI), position= position_dodge(width = 1), size = 0.5, width = 0.2)+
  facet_wrap(~ Order) +
  xlab("") +
  scale_y_continuous(limits = c(-100, 120), breaks = scales::pretty_breaks(n = 10)) +
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


ggsave(filename = paste0(outDir, "FIGURE_4_Comparison_LUSTA_Abun.pdf"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "FIGURE_4_Comparison_LUSTA_Abun.jpeg"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 600)

### now richness ###
plot_data <- final_res[final_res$Metric == "Species richness" | final_res$Metric == "Species Richness",]

# create point and error bar plot
ggplot(data = plot_data, aes(col = LUI, group = STA)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  geom_point(aes(x = LUI, y = Median, shape = STA), size = 1.5, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = LUI, ymin = Lower_CI, ymax = Upper_CI), position= position_dodge(width = 1), size = 0.5, width = 0.2)+
  facet_wrap(~ Order) +
  xlab("") +
  scale_y_continuous(limits = c(-100, 270), breaks = scales::pretty_breaks(n = 10)) +
  ylab("Percentage change in species richness (%)") +
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


ggsave(filename = paste0(outDir, "FIGURE_3_Comparison_LUSTA_Rich.pdf"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "FIGURE_3_Comparison_LUSTA_Rich.jpeg"), plot = last_plot(), width = 150, height = 120, units = "mm", dpi = 600)

