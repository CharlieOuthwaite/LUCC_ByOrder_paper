############################################################
#                                                          #
#                   Land use only models                   #
#                                                          #
############################################################

# This script takes the processed PREDICTS data and runs models of the 
# impact of land use and Order on insect biodiversity. 


# ensure working directory is clear
rm(list = ls())

# set up directories
inDir <- "1_CheckPrepareData/"
outDir <- "2_RunSimpleLUIModel/"
predsDir <- "7_Predictions/"
dataDir <- "Data/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(predsDir)) dir.create(predsDir)


# load libraries
packages <- c("StatisticalModels", "ggplot2", "cowplot", "sjPlot","dplyr", "gt", "webshot")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")


##%######################################################%##
#                                                          #
####                 1. Organise data                   ####
#                                                          #
##%######################################################%##


# read in the Site data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds")) # 7568 rows

## Species Richness Model ##

# remove NAs in the specified columns
model_data_sr <- na.omit(sites[ , c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order', 'Latitude', 'Longitude')]) # 7568 rows

# order data
model_data_sr$LUI <- factor(model_data_sr$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr$Order <- factor(model_data_sr$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"))

#save model_data_sr
saveRDS(object = model_data_sr, file = paste0(outDir,"model_data_sr.rds"))

# summaries
length(unique(model_data_sr$SS)) # 249
length(unique(model_data_sr$SSBS)) # 5642

table(model_data_sr$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       1943         581         579        3188        1277 

# look at the spread of land use/use intensity categories
table(model_data_sr$LUI)
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               1895                 1721                 1529                 2423 

table(model_data_sr$Order, model_data_sr$LUI)
#             Primary vegetation Secondary vegetation Agriculture_Low Agriculture_High
# Coleoptera                 832                  311             371              429
# Diptera                     95                  164             107              215
# Hemiptera                   94                   98             216              171
# Hymenoptera                564                  656             559             1409
# Lepidoptera                310                  492             276              199


##%######################################################%##
#                                                          #
####            2. Species Richness models              ####
#                                                          #
##%######################################################%##


# Run species richness models using GLMER function from StatisticalModels
# land use alone
sm3 <- GLMER(modelData = model_data_sr, 
             responseVar = "Species_richness", 
             fitFamily = "poisson",
             fixedStruct = "LUI", 
             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
             REML = FALSE, 
             saveVars = c('Latitude', 'Longitude'))

# save model without interaction with order
save(sm3, file = paste0(outDir, "Richness_landuse_model_noOrder.rdata"))


# check the R2 values
R2GLMER(sm3$model)
# $conditional
# [1] 0.6697481
# 
# $marginal
# [1] 0.01046651

# interaction order and LUI
sm3.3 <- GLMER(modelData = model_data_sr,
               responseVar = "Species_richness",
               fitFamily = "poisson",
               fixedStruct = "Order*LUI",
               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
               REML = FALSE, 
               saveVars = c('Latitude', 'Longitude'))

# save the model
save(sm3.3, file = paste0(outDir, "Richness_landuse_model.rdata"))

# check the R2 values
R2GLMER(sm3.3$model) # with order level interaction
# $conditional
# [1] 0.7005782
# 
# $marginal
# [1] 0.05033525

# # take a look at the AICs
AIC_sr <- AIC(sm3$model, sm3.3$model)
#             df      AIC
# sm3$model    7 46057.05
# sm3.3$model 23 43021.08

write.csv(AIC_sr, file = paste0(outDir,"AIC_sr.csv"))




##%######################################################%##
#                                                          #
####                3. Abundance Models                 ####
#                                                          #
##%######################################################%##


# remove NAs in the specified columns
model_data_ab <- na.omit(sites[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order', 'Latitude', 'Longitude')])
# 7176 rows

# order data
model_data_ab$LUI <- factor(model_data_ab$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab$Order <- factor(model_data_ab$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"))

#save model_data_ab
saveRDS(object = model_data_ab ,file = paste0(outDir,"model_data_ab.rds"))

# summaries
length(unique(model_data_ab$SS)) # 231 Studies
length(unique(model_data_ab$SSBS))# 5334

table(model_data_ab$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       1923         561         559        3006        1127 

# look at the spread of land use/use intensity categories
table(model_data_ab$LUI)
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               1803                 1583                 1510                 2280 

table(model_data_ab$Order, model_data_ab$LUI)
#             Primary vegetation Secondary vegetation Agriculture_Low Agriculture_High
# Coleoptera                 832                  311             371              409
# Diptera                     95                  164             107              195
# Hemiptera                   94                   98             216              151
# Hymenoptera                506                  594             543             1363
# Lepidoptera                276                  416             273              162



# Run abundance models using 'GLMER' function from StatisticalModels

am3 <- GLMER(modelData = model_data_ab,
             responseVar = "LogAbund",
             fitFamily = "gaussian",
             fixedStruct = "LUI",
             randomStruct = "(1|SS)+(1|SSB)",
             REML = FALSE, 
             saveVars = c('Latitude', 'Longitude'))

# save the model without an interaction with Order
save(am3, file = paste0(outDir, "Abundance_landuse_model_noOrder.rdata"))



# check the R2 values
R2GLMER(am3$model)
# $conditional
# [1] 0.3501162
# 
# $marginal
# [1] 0.03449465

# run the model with the interaction with Order
am3.3 <- GLMER(modelData = model_data_ab,
               responseVar = "LogAbund",
               fitFamily = "gaussian",
               fixedStruct = "Order*LUI",
               randomStruct = "(1|SS)+(1|SSB)",
               REML = FALSE, 
               saveVars = c('Latitude', 'Longitude'))

# save model ouputs
save(am3.3, file = paste0(outDir, "Abundance_landuse_model.rdata"))


# abundance models - no warnings #

# check the R2 values
R2GLMER(am3.3$model)
# $conditional
# [1] 0.3796676
# 
# $marginal
# [1] 0.05109883

# # take a look at the AICs
AIC_ab <- AIC(am3$model, am3.3$model)
#             df      AIC
# am3$model    7 22309.66
# am3.3$model 23 22195.28

write.csv(AIC_ab, file = paste0(outDir,"AIC_ab.csv"))



#### 4. Table of AICs ####

# species richness and abundance together
selection_table <- data.frame("Response" = c(rep("Species richness", 2),
                                             rep("Abundance", 2)),
                              "Model" = c("Species richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Abundance ~ LUI + (1|SS) + (1|SSB)",
                                          "Abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                              "AIC" = c(AIC(sm3$model), AIC(sm3.3$model), AIC(am3$model), AIC(am3.3$model)),
                              "Conditional_Rsquared" = c(R2GLMER(sm3$model)[[1]], R2GLMER(sm3.3$model)[[1]], R2GLMER(am3$model)[[1]], R2GLMER(am3.3$model)[[1]]),
                              "Marginal_Rsquared" = c(R2GLMER(sm3$model)[[2]], R2GLMER(sm3.3$model)[[2]], R2GLMER(am3$model)[[2]], R2GLMER(am3.3$model)[[2]])) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  dplyr::select(Model, AIC, deltaAIC, Conditional_Rsquared, Marginal_Rsquared) %>%
  gt(rowname_col = "Model") %>%
  tab_row_group(
    label = "Species Richness",
    rows = starts_with("Species richness")
  ) %>%
  tab_row_group(
    label = "Abundance",
    rows = starts_with("Abundance")
  )%>% 
  cols_align(
    align = "center",
    columns = c(Model, AIC, deltaAIC)
  )%>%
  tab_stubhead(label = "Models")

# save it
gtsave(selection_table, "TableSX_RSquaredAIC.png", path = outDir)


# save model output tables for use in supplementary information 
# use function from sjPlot library to save neat versions of model output table
tab_model(am3.3$model, transform = NULL, file = paste0(outDir,"SupptabX_Output_table_abund.html"))
summary(am3.3$model) # check the table against the outputs
R2GLMER(am3.3$model) # check the R2 values 
# $conditional
# [1] 0.3796676
# 
# $marginal
# [1] 0.05109883

tab_model(sm3.3$model, transform = NULL, file = paste0(outDir,"SupptabX_Output_table_rich.html"))
summary(sm3.3$model) # check the table against the outputs
R2GLMER(sm3.3$model) # check the R2 values 
$conditional
# [1] 0.7005782
# 
# $marginal
# [1] 0.05033525


# save model output tables for use in supplementary information 
# use function from sjPlot library to save neat versions of model output table
tab_model(am3$model, transform = NULL, file = paste0(outDir,"Output_table_abund_noOrder.html"))
summary(am3$model) # check the table against the outputs
R2GLMER(am3$model) # check the R2 values 
# $conditional
# [1] 0.3501162
# 
# $marginal
# [1] 0.03449465

tab_model(sm3$model, transform = NULL, file = paste0(outDir,"Output_table_rich_noOrder.html"))
summary(sm3$model) # check the table against the outputs
R2GLMER(sm3$model) # check the R2 values 
# $conditional
# [1] 0.6697481
# 
# $marginal
# [1] 0.01046651


##%######################################################%##
#                                                          #
####           Richness and abundance plots             ####
#                          Global                          #
#                                                          #
##%######################################################%##

# These plots are not presented in the paper

####  5. Species Richness Plot ####
richness_metric <- predict_effects(iterations = 1000,
                                   model = sm3.3$model,
                                   model_data = model_data_sr,
                                   response_variable = "Species_richness",
                                   fixed_number = 2,
                                   fixed_column = c("Order", "LUI"),
                                   factor_number_1 = 5,
                                   factor_number_2 = 4,
                                   neg_binom = FALSE)

# rename prediction data frame and drop "Species_richness" column
result.sr <- fin_conf
result.sr <- dplyr::select(result.sr,-c(Species_richness))

richness_metric

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_sr$LUI), 5)[1:20]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric)

####  6. Abundance Plot ####

abundance_metric <- predict_effects(iterations = 1000,
                                    model = am3.3$model,
                                    model_data = model_data_ab,
                                    response_variable = "LogAbund",
                                    fixed_number = 2,
                                    fixed_column = c("Order", "LUI"),
                                    factor_number_1 = 5,
                                    factor_number_2 = 4,
                                    neg_binom = FALSE)

# rename prediction data frame and drop "Abundance" column
result.ab <- fin_conf
result.ab <- dplyr::select(result.ab,-c(LogAbund))


abundance_metric

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab$LUI), 5)[1:20]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric)


####  7. Plot Together  ####

richness <- richness_metric +
  labs(y ="Species richness diff. (%)", x = NULL) +
  guides(scale = "none") +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50), limits = c(-100, 50)) +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        panel.border = element_rect(linewidth = 0.2))+
  ggtitle("a")


abundance <- abundance_metric +
  ylab("Total abundance diff. (%)") +
  xlab("") +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50), limits = c(-100, 50)) +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
        axis.text.y = element_text(size = 7),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        panel.border = element_rect(linewidth = 0.2))+
  ggtitle("b")


# plot together
plot_figure <- cowplot::plot_grid(richness, abundance, ncol = 1, rel_heights = c(0.8,1))

# save plot (pdf)
ggsave(filename = paste0(outDir, "SimpleLUI.pdf"), plot = last_plot(), width = 200, height = 250, units = "mm", dpi = 300)


#### 8. Table of predicted values ####

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)
all_res$measure <- c(rep("ab", 20), rep("sr", 20))

# save as .csv
write.csv(all_res, file = paste0(outDir,"LUI_predictions.csv"))


############################################################
#                                                          #
#             Plots for models excluding Order             #
#                                                          #
############################################################

# create dataframe for values to predict response to
nd <- data.frame(LUI = factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), 
                              levels = levels(sm3$data$LUI)),
                 Species_richness = 0,
                 LogAbund = 0)

## species richness predictions ##
s.preds <- StatisticalModels::PredictGLMERRandIter(model = sm3$model, data = nd)

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
a.preds <- sweep(x = a.preds, MARGIN = 2, STATS = a.preds[1,], FUN = '/')

# get quantiles
a.preds.median <- ((apply(X = a.preds,MARGIN = 1,FUN = median))*100)-100
a.preds.upper <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
a.preds.lower <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


# combine data into one table for plotting
abun_res <- as.data.frame(cbind(a.preds.median, a.preds.lower, a.preds.upper))
rich_res <- as.data.frame(cbind(s.preds.median, s.preds.lower, s.preds.upper))
colnames(abun_res) <- c("median", "lower", "upper")
colnames(rich_res) <- c("median", "lower", "upper")
abun_res$metric <- "abun"
rich_res$metric <- "rich"
abun_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))
rich_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))

abun_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA
rich_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA


# point points and error bars

p1 <- ggplot(data = abun_res) +
  geom_point(aes(x = LU, y = median, col = LU), size = 1.2) +
  geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  xlab("") +
  ylab("Change in total abundance (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
  ylim(c(-100, 0)) +
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
        axis.ticks = element_line(size = 0.2),
        axis.line = element_line(size = 0.2)) +
  ggtitle("b")


p2 <- ggplot(data = rich_res) +
  geom_point(aes(x = LU, y = median, col = LU), size = 1.2) +
  geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  xlab("") +
  ylab("Change in species richness (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
  ylim(c(-100, 0)) +
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
        axis.ticks = element_line(size = 0.2),
        axis.line = element_line(size = 0.2)) +
  ggtitle("a")


p3 <- cowplot::plot_grid(p2, p1)

# save
ggsave(filename = paste0(outDir, "FIG_No_Order_LandUseOnly.pdf"), plot = last_plot(), width = 160, height = 90, units = "mm", dpi = 300)







