############################################################
#                                                          #
#            Models without rescaling abundance            #
#                                                          #
############################################################

# In this script, the models for land use and land use/climate impacts
# are rerun for total abundance, however abundance values are not rescaled
# within study/order. 

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


##%######################################################%##
#                                                          #
####                Land use only models                ####
#                                                          #
##%######################################################%##

# read in the Site data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds")) # 8884 rows

# remove NAs in the specified columns
model_data_ab <- na.omit(sites[,c('LogAbund_noRS','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order', 'Latitude', 'Longitude')])
# 8496 rows

# order data
model_data_ab$LUI <- factor(model_data_ab$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab$Order <- factor(model_data_ab$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"))

# Run the abundance model using the log abundance data that has not been rescaled 
am3.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund_noRS",fitFamily = "gaussian",
               fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE, saveVars = c('Latitude', 'Longitude'))

# save model output
save(am3.3, file = paste0(outDir, "Abundance_landuse_model_noRS.rdata"))


# check the R2 values
R2GLMER(am3.3$model)
# $conditional
# [1] 0.6026711
# 
# $marginal
# [1] 0.07881955

library(sjPlot)
# save model output tables for use in supplementary information 
# use function from sjPlot library to save neat versions of model output table
tab_model(am3.3$model, transform = NULL, file = paste0(outDir,"SupptabX_Output_table_abund_noRS.html"))
summary(am3.3$model) # check the table against the outputs
R2GLMER(am3.3$model) # check the R2 values 
# $conditional
# [1] 0.602671
# 
# $marginal
# [1] 0.07881956



#### Abundance Plot ####

abundance_metric <- predict_effects(iterations = 1000,
                                    model = am3.3$model,
                                    model_data = model_data_ab,
                                    response_variable = "LogAbund_noRS",
                                    fixed_number = 2,
                                    fixed_column = c("Order", "LUI"),
                                    factor_number_1 = 5,
                                    factor_number_2 = 4,
                                    neg_binom = FALSE)

# rename prediction data frame and drop "Abundance" column
result.ab <- fin_conf
result.ab <- dplyr::select(result.ab,-c(LogAbund_noRS))

# adjust plot
abundance <- abundance_metric +
  labs(y ="Total abundance diff. (%)", x = "Order") +
  scale_y_continuous(limits = c(-100, 140)) + 
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "right")


# save plot (pdf)
ggsave(filename = paste0(outDir, "SimpleLUI_ab_noRS.pdf"), plot = last_plot(), width = 150, height = 100, units = "mm", dpi = 300)

# function to extract prediction data
model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab$LUI), 5)[1:20]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
ab_pred <- model_data(abundance_metric)

# save prediction data
write.csv(ab_pred, file = paste0(outDir, "/Predictions_simplemod_ab_noRS.csv"), row.names = F)



############################################################
#                                                          #
####               LUI and STA models                   ####
#                                                          #
############################################################


# load in the data with climate variable
predictsSites <- readRDS(file = paste0(dataDir,"PREDICTSSitesClimate_Data.rds"))

# i. Abundance, mean anomaly including interaction
# subset to those sites with abundance data
model_data <- predictsSites[!is.na(predictsSites$LogAbund_noRS), ]

MeanAnomalyModelAbund_noresc <- GLMER(modelData = model_data, 
                                      responseVar = "LogAbund_noRS",
                                      fitFamily = "gaussian",
                                      fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                                      randomStruct = "(1|SS)+(1|SSB)",
                                      saveVars = c("SSBS", 'Latitude', 'Longitude'))

# get summary
summary(MeanAnomalyModelAbund_noresc$model)


# save the model output
save(MeanAnomalyModelAbund_noresc, file = paste0(outDir, "LUICCmodel_abun_norescaling.rdata"))

tab_model(MeanAnomalyModelAbund_noresc$model, transform = NULL, file = paste0(outDir,"SupptabX_Output_table_abund_CC_noRS.html"))
summary(MeanAnomalyModelAbund_noresc$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_noresc$model) # check the R2 values 
# $conditional
# [1] 0.6593063
# 
# $marginal
# [1] 0.1331049


##### predictions for LUI CC model #####


# what is the rescaled value of STA of 1
BackTransformCentreredPredictor(transformedX = 0.999, originalX = predictsSites$StdTmeanAnomaly) # 0.999 gives about 1 

# what is the rescaled value of STA of 0
BackTransformCentreredPredictor(transformedX = -1.39, originalX = predictsSites$StdTmeanAnomaly) # -1.39 gives about 0 

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
nd$LogAbund_noRS <- 0

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


final_res <- result.ab

final_res$Order <- factor(final_res$Order, levels = c("All insects", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"))

final_res$STA <- c(0, 1)
final_res$STA <- factor(final_res$STA, levels = c("0", "1"))

final_res$LUI <- factor(final_res$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

names(final_res)[1:3] <- c("Median",  "Upper_CI", "Lower_CI")


# create point and error bar plot
ggplot(data = final_res, aes(col = LUI, group = STA)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  geom_point(aes(x = LUI, y = Median, shape = STA), size = 1.5, position= position_dodge(width = 1)) + 
  geom_errorbar(aes(x = LUI, ymin = Lower_CI, ymax = Upper_CI), position= position_dodge(width = 1), size = 0.5, width = 0.2)+
  facet_wrap(~ Order, scales = "free") +
  xlab("") +
  #scale_y_continuous(limits = c(-100, 150), breaks = scales::pretty_breaks(n = 10)) +
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

# save the plot
ggsave(filename = paste0(outDir, "SuppFIGX_LUSTA_Abun_noRS.pdf"), plot = last_plot(), width = 120, height = 150, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "SuppFIGX_LUSTA_Abun_noRS.jpeg"), plot = last_plot(), width = 120, height = 150, units = "mm", dpi = 600)
