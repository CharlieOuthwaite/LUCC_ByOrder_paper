##%######################################################%##
#                                                          #
####               Prepare PREDICTS data                ####
#                                                          #
##%######################################################%##

# This script takes the 2016 release of the PREDICTS database, subsets it to insects, 
# and organises the data by Order. 

# ensure working directory is clear
rm(list = ls())

# set up directories
dataDir <- "Data/"
outDir <- "1_CheckPrepareData/"
if(!dir.exists(outDir)) dir.create(outDir)
plotDir <- "1_CheckPrepareData/Plots/"
if(!dir.exists(plotDir)) dir.create(plotDir)

# Load required libraries
packages <- c("predictsFunctions","patchwork", "dplyr", "ggplot2", "yarg", "lme4", "gt",
              "broom.mixed", "MASS","webshot")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

# source in additional functions
source("0_Functions.R")

#### 1. Load and Organise data ####

# Set the path to copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)
# 3250404 obs. of 67 variables

# Select only data for insects
predicts <- predicts[(predicts$Class=="Insecta"),]
# 935078 obs. of 67 variables

# Correct effort-sensitive abundance measures (assumes linear relationship between 
# effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
# Correcting 0 missing sampling effort values
# Re-scaling sampling effort
# Correcting 870378 values for sensitivity to sampling effort 

# check diversity metrics
table(predicts$Diversity_metric)

# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those entries that are the problem
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]
# 934845 obs. of 67 variables

# Merge sites that have the same coordinates, from the same study and same taxonomic 
# family (e.g. multiple traps on a single transect)
predicts <- predictsFunctions::MergeSites(diversity = predicts)
# 826292 obs. of 67 variables

# remove entries without Order
predicts <- droplevels(predicts[(predicts$Order!=""),])
# 826016 obs. of 67 variables

# check
table(predicts$Order) 

# keep top five orders 
predicts <- predicts %>% filter(Order %in% c("Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera", "Hemiptera")) %>% droplevels()
# 810399 obs. of 67  variables

# convert Order to a "factor"
predicts$Order <- as.factor(predicts$Order)

#### 2. Calculate site metrics and prepare dataset ####

# Split predicts into separate data frames according to insect Order 

# use split function to split the predicts data frame into 5 data frames (1/Order)
OrderName <- paste0("",predicts$Order)

by_Order <- split(predicts,OrderName)

# extract data frames from list into global environment
list2env(by_Order,globalenv())

# drop levels and calculate site level metrics
Coleoptera <- droplevels(Coleoptera)
Coleoptera <- SiteMetrics(diversity = Coleoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))

Diptera <- droplevels(Diptera)
Diptera <- SiteMetrics(diversity = Diptera,
                       extra.cols = c("Predominant_land_use",
                                      "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))

Hemiptera <- droplevels(Hemiptera)
Hemiptera <- SiteMetrics(diversity = Hemiptera,
                         extra.cols = c("Predominant_land_use",
                                        "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))

Hymenoptera <- droplevels(Hymenoptera)
Hymenoptera <- SiteMetrics(diversity = Hymenoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))

Lepidoptera <- droplevels(Lepidoptera)
Lepidoptera <- SiteMetrics(diversity = Lepidoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))

# merge all sites_Order data frames into one called "sites"
sites <- rbind(Coleoptera, Diptera, Hemiptera, Hymenoptera, Lepidoptera)
# 10746 rows

# First, we will rearrange the land-use classification 
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

# 10746 obs. of 29 variables

# remove the urban sites and sites that are NA in LUI
sites <- sites[!sites$LUI == "Urban", ]
sites <- sites[!is.na(sites$LUI), ]
# 8890 obs. of 29 variables

sites <- droplevels(sites)

# rescale abundance values
sites <- RescaleAbundance(sites)

# add another column without rescaling for testing later
# transform abundance values 
sites$LogAbund_noRS <- log(sites$Total_abundance+0.01)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]
# 8884 obs. of 30 variables

# create a new variable designating sites as Tropical or Non-tropical
# assign new variable for tropical/temperate, convert to factor, and filter out NA
sites$Realm <- ifelse(sites$Latitude >= -23.5 & sites$Latitude <= 23.5, "Tropical", "NonTropical")
sites$Realm <- factor(sites$Realm, levels = c("NonTropical", "Tropical"))
sites <- sites %>%
  filter(!is.na(Realm))

# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds")) 

##%######################################################%##
#                                                          #
####             3. Summarise data by order             ####
#                                                          #
##%######################################################%##


# n unique sites across the dataset
length(unique(sites$SSBS)) # 6014
# n unique studies across the dataset
length(unique(sites$SS)) # 254


# This creates a summary table splitting by Order and land use to show n studies and n sites
sites_summary <- sites %>%
  group_by(Order) %>% 
  #mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
  summarise(Sites = length(SSBS),
            Studies = n_distinct(SS),
            Species = n_distinct(Best_guess_binomial),
            sites_1 = length(LUI[LUI == "Primary vegetation"]) ,
            studies_1 = n_distinct(SS[LUI == "Primary vegetation"]),
            sites_2 = length(LUI[LUI == "Secondary vegetation"]),
            studies_2 = n_distinct(SS[LUI == "Secondary vegetation"]),
            sites_low = length(LUI[LUI == "Agriculture_Low"]),
            studies_low = n_distinct(SS[LUI == "Agriculture_Low"]),
            sites_high = length(LUI[LUI == "Agriculture_High"]),
            studies_high = n_distinct(SS[LUI == "Agriculture_High"]),
            Abundance = sum(!is.na(LogAbund)) ,
            SpeciesRichness = sum(Species_richness>=0)) %>%
  ungroup() %>%  
  arrange(desc(Sites))%>%
  gt() %>%
  # tab_header(title = "Data spread across land-uses"
  # ) %>%
  tab_spanner(
    label = "Total",
    columns = c(Sites,Studies,Species)
  )  %>%
  tab_spanner(
    label = "Primary Vegetation",
    columns = c(sites_1,studies_1)
  ) %>%
  tab_spanner(
    label = "Secondary vegetation",
    columns = c(sites_2,studies_2)
  ) %>%
  tab_spanner(
    label = "Low-intensity agriculture",
    columns = c(sites_low,studies_low)
  ) %>%
  tab_spanner(
    label = "High-intensity agriculture",
    columns = c(sites_high,studies_high)
  )  %>%
  tab_spanner(
    label = "Diversity Metric",
    columns = c(Abundance,SpeciesRichness)
    )  %>% 
  cols_align(
    align = "center",
    columns = c(Order,Sites,Studies,Species,sites_1,studies_1,sites_2,studies_2,sites_low,studies_low,sites_high,studies_high,Abundance,SpeciesRichness)
  ) %>%
  cols_label(
    Order = "Order",
    Sites = "Sites",
    Studies = "Studies",
    Species = "Species",
    sites_1 = "Sites",
    studies_1 = "Studies",
    sites_2 = "Sites",
    studies_2 = "Studies",
    sites_low = "Sites",
    studies_low = "Studies",
    sites_high = "Sites",
    studies_high = "Studies",
    Abundance = "Abundance",
    SpeciesRichness = "Species richness"
  )

# save table
gtsave(sites_summary,"sites_summary_SuppTab1.png",path = outDir)

# save as csv
write.csv(sites_summary,file = paste0(outDir, "sites_summary_SuppTab1.csv"))

# assess the record level dataset
# subset to those sites and orders selected
pred_sub <- predicts[predicts$SSBS %in% sites$SSBS & predicts$Order %in% sites$Order, ]
pred_sub <- droplevels(pred_sub)
nrow(pred_sub) # 741933 records
table(pred_sub$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#     377237       24149       45217      144474      150856

table(pred_sub$Predominant_land_use)
length(unique(pred_sub$SS)) # 254
length(unique(pred_sub$SSBS)) # 6014
length(unique(pred_sub$Taxon_name_entered)) # 17218 "species"

# determine unique taxa counts
species <- unique(pred_sub[,c('Order',"Family", 'Taxon_name_entered')])

# order level counts of taxon names supplied
order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       6037        1477        1302        4528        3874


# Order level counts of n families
family.counts <- tapply(X = species$Family,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#        103          76          68          54          60


############################################################
#                                                          #
####    Figure 1: map and bar chart of sites/orders     ####
#                                                          #
############################################################

# load the prepared dataset if not already loaded
sites <- readRDS(file = paste0(outDir,"PREDICTSSiteData.rds")) # 8884 rows

# load map of country borders using rgdal
map.borders1 <- readOGR(dsn = paste0(dataDir,"World Map/TM_WORLD_BORDERS-0.3.shp"), #Provide the directory of shp file 
                        #layer = "TM_WORLD_BORDERS-0.3", #Provide the name of the shp file without extension (.shp)
                        verbose = FALSE)

# convert to dataframe where the id column will be SUBREGION from the @data
shp_df <- broom::tidy(map.borders1, region = "SUBREGION")

## 1. plot map of site distribution: scatterplot ##

# plot the raster in ggplot

# map of sites, no colour
p_map1.1 <- ggplot() +
  geom_polygon(data = shp_df, aes(x = long, y = lat, group = group), alpha = 0.5, col = "black", fill = c("grey")) +
  coord_fixed() +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, col = Order), shape = 21, size=2, col = c("#A4D3EE"), fill = c("#104E8B")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")


## Organise sites and order by UN subregion
## UN Subregion ##

# need to make a new dataframe with variables:
# UN_subregion
# Latitude of centre of subregion
# Longitude of centre of subregion
# Order
# Sum of records for each order in that subregion

mapdf <- sites %>% dplyr::select(UN_subregion,Order,Longitude,Latitude) %>% # filter for rows I want
  group_by(UN_subregion) %>% # group by subregion
  summarise(Coleoptera = length(Order[Order == "Coleoptera"]), # count for each subregion, the number of records per order
            Diptera = length(Order[Order == "Diptera"]),
            Hemiptera = length(Order[Order == "Hemiptera"]),
            Hymenoptera = length(Order[Order == "Hymenoptera"]),
            Lepidoptera = length(Order[Order == "Lepidoptera"]),
            radius = 25*sqrt(length(UN_subregion)/nrow(sites))) # calculate the radius for each subregion


## create a barplot showing n sites in each subregion ##

# need to reorganize dataframe for barplot
mapdf2 <- mapdf %>% dplyr::select(UN_subregion, Coleoptera, Diptera, Hemiptera, Hymenoptera, Lepidoptera)
mapdf2 <- cbind(mapdf2[1], stack(mapdf2[2:3:4:5:6]))
mapdf2 <- mapdf2 %>% rename_at('ind', ~'Order')
mapdf2 <- dplyr::select(mapdf2,UN_subregion, Order, values)

mapdf2$UN_subregion <- factor(mapdf2$UN_subregion)
mapdf2$Order <- factor(mapdf2$Order) 

p_bar <- ggplot(data=mapdf2) +
  scale_fill_manual(values=c("#0F6B99FF","#FFD8B2FF","#A3CC51FF","#8F7EE5FF","#990F0FFF"))+
  geom_bar(stat="identity",aes(x=UN_subregion, y=values, fill=Order), alpha = 0.7)+
  labs(x = "UN subregion", y = "Number of sites") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(0, 1600), breaks = c(0, 250, 500, 750, 1000, 1250, 1500))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 50, size = 10, vjust = 0.5, hjust = 0.7),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank())


# put map and bar chart together
scatter_bar_mercator <- cowplot::plot_grid(p_map1.1,p_bar, ncol = 1, labels = c("(a)", "(b)"))

# save plot (jpeg)
ggsave("FIGURE_1_scatter_bar_mercator_UN.jpeg", device ="jpeg", path = plotDir, width=250, height=300, units="mm", dpi = 350)




