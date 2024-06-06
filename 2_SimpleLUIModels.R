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
packages_model <- c("StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "lme4", "gt", "broom.mixed", "MASS","webshot")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")


##%######################################################%##
#                                                          #
####                 1. Organise data                   ####
#                                                          #
##%######################################################%##


# read in the Site data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds")) # 8884 rows

## Species Richness Model ##

# remove NAs in the specified columns
model_data_sr <- na.omit(sites[,c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order', 'Latitude', 'Longitude')]) # 8884 rows

# order data
model_data_sr$LUI <- factor(model_data_sr$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr$Order <- factor(model_data_sr$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"))

#save model_data_sr
saveRDS(object = model_data_sr ,file = paste0(outDir,"model_data_sr.rds"))

# summaries
length(unique(model_data_sr$SS)) # 254
length(unique(model_data_sr$SSBS)) # 6014

table(model_data_sr$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       2315         721         871        3328        1649

# look at the spread of land use/use intensity categories
table(model_data_sr$LUI)
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               2093                 2392                 1886                 2513  

table(model_data_sr$Order, model_data_sr$LUI)
#             Primary vegetation Secondary vegetation Agriculture_Low Agriculture_High
# Coleoptera                 898                  504             466              447
# Diptera                     95                  250             143              233
# Hemiptera                  160                  211             311              189
# Hymenoptera                564                  742             595             1427
# Lepidoptera                376                  685             371              217


##%######################################################%##
#                                                          #
####            2. Species Richness models              ####
#                                                          #
##%######################################################%##


# Run species richness models using GLMER function from StatisticalModels
# land use alone
sm3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE, saveVars = c('Latitude', 'Longitude'))

# save model without interaction with order
save(sm3, file = paste0(outDir, "Richness_landuse_model_noOrder.rdata"))


# check the R2 values
R2GLMER(sm3$model)
# $conditional
# [1] 0.6700182
# 
# $marginal
# [1] 0.007716681

# interaction order and LUI
sm3.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE, saveVars = c('Latitude', 'Longitude'))

save(sm3.3, file = paste0(outDir, "Richness_landuse_model.rdata"))

# check the R2 values
R2GLMER(sm3.3$model) # with order level interaction
# $conditional
# [1] 0.7106983
# 
# $marginal
# [1] 0.05387397


# richness models - no warnings

# # # Run richness models using 'glmr' function from lme4
# # # need these to run allFit()
# # # run allFit() to test the models
# 
# g_sm3 <- glmer(Species_richness~LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
# g_sm3.3 <- glmer(Species_richness~Order*LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
# 
# AIC(g_sm3, g_sm3.3)
# #         df      AIC
# # g_sm3    7 64132.71
# # g_sm3.3 23 59449.28
# 
# # Testing richness models with allFit()
# 
# g_sm3_all <- allFit(g_sm3)
# g_sm3.3_all <- allFit(g_sm3.3)



##%######################################################%##
#                                                          #
####                3. Abundance Models                 ####
#                                                          #
##%######################################################%##


# remove NAs in the specified columns
model_data_ab <- na.omit(sites[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order', 'Latitude', 'Longitude')])
# 8492 rows

# order data
model_data_ab$LUI <- factor(model_data_ab$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab$Order <- factor(model_data_ab$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"))

#save model_data_ab
saveRDS(object = model_data_ab ,file = paste0(outDir,"model_data_ab.rds"))

# summaries
length(unique(model_data_ab$SS)) # 236 Studies
length(unique(model_data_ab$SSBS))# 5706

table(model_data_ab$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       2295         701         851        3146        1499

# look at the spread of land use/use intensity categories
table(model_data_ab$LUI)
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               2001                 2254                  1867                2370

table(model_data_ab$Order, model_data_ab$LUI)
#             Primary vegetation Secondary vegetation Agriculture_Low Agriculture_High
# Coleoptera                 898                  504             466              427
# Diptera                     95                  250             143              213
# Hemiptera                  160                  211             311              169
# Hymenoptera                506                  680             579             1381
# Lepidoptera                342                  609             368              180

# Run abundance models using 'GLMER' function from StatisticalModels

# am0 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#              fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am1 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#              fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#              fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE, saveVars = c('Latitude', 'Longitude'))

# save the model without an interaction with Order
save(am3, file = paste0(outDir, "Abundance_landuse_model_noOrder.rdata"))



# check the R2 values
R2GLMER(am3$model)
# $conditional
# [1] 0.3364167
# 
# $marginal
# [1] 0.02482347


# am4 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#              fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am0.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#                fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am1.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#                fixedStruct = "Order+LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am2.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#                fixedStruct = "Order+LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am3.2<- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#               fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am4.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#                fixedStruct = "Order+LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am1.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#                fixedStruct = "Order*LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# 
# am2.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#                fixedStruct = "Order*LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE, saveVars = c('Latitude', 'Longitude'))

# am4.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
#                fixedStruct = "Order*LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
# warning

save(am3.3, file = paste0(outDir, "Abundance_landuse_model.rdata"))


# abundance models - no warnings

# check the R2 values
R2GLMER(am3.3$model)
# $conditional
# [1] 0.3775368
# 
# $marginal
# [1] 0.04830284

# # take a look at the AICs
# AIC_am <-print(AIC(am0$model,am1$model,am2$model,am3$model,am4$model,
#                    am0.2$model,am1.2$model,am2.2$model,am3.2$model,am4.2$model,
#                    am1.3$model,am2.3$model,am3.3$model,am4.3$model))

AIC(am3$model, am3.3$model)

write.csv(AIC_am, file = paste0(outDir,"AIC_am.csv"))

# save 
# saveRDS(object = sm0 ,file = paste0(outDir,"sm0.rd"))
saveRDS(object = sm3 ,file = paste0(outDir,"sm3.rds"))
# saveRDS(object = sm0.2 ,file = paste0(outDir,"sm0.2.rds"))
# saveRDS(object = sm3.2 ,file = paste0(outDir,"sm3.2.rds"))
saveRDS(object = sm3.3 ,file = paste0(outDir,"sm3.3.rds"))
# saveRDS(object = am0 ,file = paste0(outDir,"am0.rds"))
saveRDS(object = am3 ,file = paste0(outDir,"am3.rds"))
# saveRDS(object = am0.2 ,file = paste0(outDir,"am0.2.rds"))
# saveRDS(object = am3.2 ,file = paste0(outDir,"am3.2.rds"))
saveRDS(object = am3.3 ,file = paste0(outDir,"am3.3.rds"))


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
# [1] 0.3775368
# 
# $marginal
# [1] 0.04830284

tab_model(sm3.3$model, transform = NULL, file = paste0(outDir,"SupptabX_Output_table_rich.html"))
summary(sm3.3$model) # check the table against the outputs
R2GLMER(sm3.3$model) # check the R2 values 
# $conditional
# [1] 0.7106989
# 
# $marginal
# [1] 0.05387385


# save model output tables for use in supplementary information 
# use function from sjPlot library to save neat versions of model output table
tab_model(am3$model, transform = NULL, file = paste0(outDir,"Output_table_abund_noOrder.html"))
summary(am3$model) # check the table against the outputs
R2GLMER(am3$model) # check the R2 values 
# $conditional
# [1] 0.3364167
# 
# $marginal
# [1] 0.02482347

tab_model(sm3$model, transform = NULL, file = paste0(outDir,"Output_table_rich_noOrder.html"))
summary(sm3$model) # check the table against the outputs
R2GLMER(sm3$model) # check the R2 values 
# $conditional
# [1] 0.6700185
# 
# $marginal
# [1] 0.007716608


##%######################################################%##
#                                                          #
####           Richness and abundance plots             ####
#                          Global                          #
#                                                          #
##%######################################################%##
# 
# # These plots are not presented in the paper
# 
# ####  5. Species Richness Plot ####
# richness_metric <- predict_effects(iterations = 1000,
#                                    model = sm3.3$model,
#                                    model_data = model_data_sr,
#                                    response_variable = "Species_richness",
#                                    fixed_number = 2,
#                                    fixed_column = c("Order", "LUI"),
#                                    factor_number_1 = 5,
#                                    factor_number_2 = 4,
#                                    neg_binom = FALSE)
# 
# # rename prediction data frame and drop "Species_richness" column
# result.sr <- fin_conf
# result.sr <- dplyr::select(result.sr,-c(Species_richness))
# 
# richness_metric
# 
# model_data <- function(model_plot){
#   ggplot_build(model_plot)$data[[2]] %>%
#     dplyr::select(y) %>%
#     cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
#     mutate(LUI = rep(levels(model_data_sr$LUI), 5)[1:20]) %>%
#     dplyr::select(LUI, y, ymin, ymax)
# }
# 
# # richness data
# model_data(richness_metric)
# 
# ####  6. Abundance Plot ####
# 
# abundance_metric <- predict_effects(iterations = 1000,
#                                     model = am3.3$model,
#                                     model_data = model_data_ab,
#                                     response_variable = "LogAbund",
#                                     fixed_number = 2,
#                                     fixed_column = c("Order", "LUI"),
#                                     factor_number_1 = 5,
#                                     factor_number_2 = 4,
#                                     neg_binom = FALSE)
# 
# # rename prediction data frame and drop "Abundance" column
# result.ab <- fin_conf
# result.ab <- dplyr::select(result.ab,-c(LogAbund))
# 
# 
# abundance_metric
# 
# model_data <- function(model_plot){
#   ggplot_build(model_plot)$data[[2]] %>%
#     dplyr::select(y) %>%
#     cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
#     mutate(LUI = rep(levels(model_data_ab$LUI), 5)[1:20]) %>%
#     dplyr::select(LUI, y, ymin, ymax)
# }
# 
# # abundance data
# model_data(abundance_metric)
# 
# 
# ####  7. Plot Together  ####
# 
# richness <- richness_metric + 
#   labs(y ="Species richness diff. (%)", x = NULL) +
#   guides(scale = "none") +
#   scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
#         axis.text.y = element_text(size = 7),
#         legend.position = "none")+
#   ggtitle("a")
# 
# 
# abundance <- abundance_metric +
#   labs(y ="Total abundance diff. (%)", x = "Order") +
#   scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
#   theme(axis.title = element_text(size = 8),
#         axis.text.x = element_text(size = 7,angle=45,margin=margin(t=20)),
#         axis.text.y = element_text(size = 7),
#         legend.position = "none")+
#   ggtitle("b")
# 
# 
# # get the legend
# legend <- get_legend(
#   abundance +
#     guides(color = guide_legend(ncol = 4)) +
#     theme(legend.position = "bottom",
#           legend.box = "vertical", 
#           legend.text = element_text(size = 7), 
#           legend.title = element_blank())
# )
# 
# # plot together
# plot_figure <- cowplot::plot_grid(richness, abundance, ncol = 1, rel_heights = c(1,1))
# legend <- cowplot::plot_grid(legend, ncol = 1, rel_heights = c(0.5,1,0.5))
# plot_figure <- cowplot::plot_grid(plot_figure, legend, nrow = 2, rel_heights = c(1,0.2))
# 
# # save plot (pdf)
# ggsave(filename = paste0(outDir, "SimpleLUI.pdf"), plot = last_plot(), width = 150, height = 200, units = "mm", dpi = 300)
# 
# # save plot (jpeg)
# ggsave("SimpleLUI.jpeg", device ="jpeg", path = outDir, width=20, height=15, units="cm", dpi = 350)
# 
# #### 8. Table of predicted values ####
# 
# # combine results into a table for saving
# all_res <- rbind(result.ab, result.sr)
# all_res$measure <- c(rep("ab", 20), rep("sr", 20))
# 
# # save as .csv
# write.csv(all_res, file = paste0(predsDir,"LUI_predictions.csv"))
# 
# # save as table
# LUI_predictions <- all_res %>% gt()
# gtsave(LUI_predictions,"LUI_predictions.png", path = predsDir)
# 
# 
# ############################################################
# #                                                          #
# #             Plots for models excluding Order             #
# #                                                          #
# ############################################################
# 
# # create dataframe for values to predict response to
# nd <- data.frame(LUI=factor(c("Primary vegetation","Secondary vegetation",
#                               "Agriculture_Low","Agriculture_High"), levels = levels(sm3$data$LUI)),
#                  Species_richness=0,
#                  LogAbund=0)
# 
# ## species richness predictions ##
# s.preds <- StatisticalModels::PredictGLMERRandIter(model = sm3$model, data = nd)
# 
# s.preds <- exp(s.preds)
# 
# # convert to percentage difference from primary vegetation
# s.preds <- sweep(x = s.preds, MARGIN = 2, STATS = s.preds[1,], FUN = '/')
# 
# # get quantiles
# s.preds.median <- ((apply(X = s.preds,MARGIN = 1,FUN = median))*100)-100
# s.preds.upper <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
# s.preds.lower <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100
# 
# 
# 
# ## abundance predictions ##
# a.preds <- PredictGLMERRandIter(model = am3$model,data = nd)
# 
# a.preds <- exp(a.preds)-0.01
# 
# # convert to percentage difference from primary vegetation
# a.preds <- sweep(x = a.preds, MARGIN = 2, STATS = a.preds[1,], FUN = '/')
# 
# # get quantiles
# a.preds.median <- ((apply(X = a.preds,MARGIN = 1,FUN = median))*100)-100
# a.preds.upper <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
# a.preds.lower <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100
# 
# 
# # combine data into one table for plotting
# abun_res <- as.data.frame(cbind(a.preds.median, a.preds.lower, a.preds.upper))
# rich_res <- as.data.frame(cbind(s.preds.median, s.preds.lower, s.preds.upper))
# colnames(abun_res) <- c("median", "lower", "upper")
# colnames(rich_res) <- c("median", "lower", "upper")
# abun_res$metric <- "abun"
# rich_res$metric <- "rich"
# abun_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))
# rich_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))
# 
# abun_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA
# rich_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA
# 
# 
# # point points and error bars
# 
# p1 <- ggplot(data = abun_res) +
#   geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
#   geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
#   geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
#   xlab("") +
#   ylab("Change in total abundance (%)") +
#   scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
#   ylim(c(-100, 0)) +
#   theme(legend.position = "none", 
#         aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text.y = element_text(size = 7),
#         axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
#         axis.title = element_text(size = 7),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(), 
#         panel.background = element_blank(), 
#         axis.ticks = element_line(size = 0.2), 
#         axis.line = element_line(size = 0.2)) +
#   ggtitle("b")
# 
# 
# p2 <- ggplot(data = rich_res) +
#   geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
#   geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
#   geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
#   xlab("") +
#   ylab("Change in species richness (%)") +
#   scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
#   ylim(c(-100, 0)) +
#   theme(legend.position = "none", 
#         aspect.ratio = 1, 
#         title = element_text(size = 8, face = "bold"),
#         axis.text.y = element_text(size = 7),
#         axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
#         axis.title = element_text(size = 7),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(), 
#         panel.background = element_blank(), 
#         axis.ticks = element_line(size = 0.2), 
#         axis.line = element_line(size = 0.2)) +
#   ggtitle("a")
# 
# 
# p3 <- cowplot::plot_grid(p2, p1)
# 
# # save
# ggsave(filename = paste0(outDir, "FIG_No_Order_LandUseOnly.pdf"), plot = last_plot(), width = 160, height = 90, units = "mm", dpi = 300)
# 






