##%######################################################%##
#                                                          #
####               Prepare PREDICTS data                ####
#                                                          #
##%######################################################%##

# This script takes the 2016 release of the PREDICTS database, subsets it to insects, 
# and organises the data by Order. 

# UPDATE: outliers that have a strong affect on the results are now removed.

# ensure working directory is clear
rm(list = ls())

# set up directories
dataDir <- "Data/"
outDir <- "1_CheckPrepareData/"
if(!dir.exists(outDir)) dir.create(outDir)
plotDir <- "1_CheckPrepareData/Plots/"
if(!dir.exists(plotDir)) dir.create(plotDir)

# Load required libraries
packages <- c("predictsFunctions", "dplyr", "ggplot2", "gt", "webshot2")
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

# what % of records are from the 5 main orders
table(droplevels(predicts$Order))
# 457507 + 31711 + 52187 + 206165 + 168650   = 916220                                    
# 916220/nrow(predicts)*100 = 97.98327 %

# how many records are at the genus, species or infraspecies level?
nrow(predicts[predicts$Rank %in% c("Infraspecies", "Species", "Genus"), ]) # 592161
# 592161/nrow(predicts)*100 = 63.32744 %

# remove the 4 studies identified as influential outliers
outliers <- c("SC1_2005__Richardson 1", "HW1_2011__Summerville 1", "SC1_2011__Meijer 1", "AD1_2008__Billeter 6")
# take a look at the outliers
outdata <- predicts[predicts$SS %in% outliers, ] # 118584 obs
length(unique(outdata$SS)) # 4
length(unique(outdata$SSBS)) # 368

# remove outliers
predicts <- predicts[!predicts$SS %in% outliers, ] # 816494 obs

# Correct effort-sensitive abundance measures (assumes linear relationship between 
# effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
# Correcting 0 missing sampling effort values
# Rescaling sampling effort
# Correcting 751794 values for sensitivity to sampling effort

# check diversity metrics
table(predicts$Diversity_metric)

# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those entries that are the problem
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]
# 816261 obs. of 67 variables

# Merge sites that have the same coordinates, from the same study and same taxonomic 
# family (e.g. multiple traps on a single transect)
predicts <- predictsFunctions::MergeSites(diversity = predicts)
# 713207 obs. of 67 variables

# remove entries without Order
predicts <- droplevels(predicts[(predicts$Order!=""),])
# 712931 obs. of 67 variables

# check
table(predicts$Order) 

# keep top five orders 
predicts <- predicts %>% filter(Order %in% c("Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera", "Hemiptera")) %>% droplevels()
# 703796 obs. of 67  variables

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
# 9707 rows

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

# 9707 obs. of 28 variables

# remove the urban sites and sites that are NA in LUI
sites <- sites[!sites$LUI == "Urban", ]
sites <- sites[!is.na(sites$LUI), ]
# 7851 obs. of 28 variables

sites <- droplevels(sites)

# rescale abundance values, this also log-transforms values, adding 0.01 to all
sites <- RescaleAbundance(sites)

# add another column without rescaling for testing later
# transform abundance values 
sites$LogAbund_noRS <- log(sites$Total_abundance+0.01)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]
# 7845 obs. of 30 variables

# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds")) 

##%######################################################%##
#                                                          #
####             3. Summarise data by order             ####
#                                                          #
##%######################################################%##


# n unique sites across the dataset
length(unique(sites$SSBS)) # 5659
# n unique studies across the dataset
length(unique(sites$SS)) # 250


# This creates a summary table splitting by Order and land use to show n studies and n sites
sites_summary <- sites %>%
  group_by(Order) %>% 
  #mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
  summarise(Sites = length(SSBS),
            Studies = n_distinct(SS),
            #Species = n_distinct(Best_guess_binomial),
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
    columns = c(Sites,Studies)
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
    columns = c(Order,Sites,Studies,sites_1,studies_1,sites_2,studies_2,sites_low,studies_low,sites_high,studies_high,Abundance,SpeciesRichness)
  ) %>%
  cols_label(
    Order = "Order",
    Sites = "Sites",
    Studies = "Studies",
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
gtsave(sites_summary,"sites_summary_SuppTab1_outrm.png", path = outDir)

# save as csv
write.csv(sites_summary,file = paste0(outDir, "sites_summary_SuppTab1_outrm.csv"))

# assess the record level dataset
# subset to those sites and orders selected
pred_sub <- predicts[predicts$SSBS %in% sites$SSBS & predicts$Order %in% sites$Order, ]
pred_sub <- droplevels(pred_sub)
nrow(pred_sub) # 635330 records
table(pred_sub$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
# 337823       18062       30690      137806      110949 

table(pred_sub$Predominant_land_use)
length(unique(pred_sub$SS)) # 250
length(unique(pred_sub$SSBS)) # 5659
length(unique(pred_sub$Taxon_name_entered)) # 16357 "species"

abun <- sites[!is.na(sites$LogAbund), ]
nrow(abun) # 7453
length(unique(abun$SS)) # 232
length(unique(abun$SSBS)) # 5351

# determine unique taxa counts
species <- unique(pred_sub[,c('Order',"Family", 'Taxon_name_entered')])

# order level counts of taxon names supplied
order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       5717        1384        1197        4470        3589  


# Order level counts of n families
family.counts <- tapply(X = species$Family,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#        101          75          64          54          59


############################################################
#                                                          #
####    Figure 1: map and bar chart of sites/orders     ####
#                                                          #
############################################################

# load the prepared dataset if not already loaded
sites <- readRDS(file = paste0(outDir,"PREDICTSSiteData.rds")) # 7848 rows

# plot the raster in ggplot
map.world <- map_data('world')


## 1. plot map of site distribution: scatterplot ##

# plot the raster in ggplot

# map of sites, no colour
p_map1.1 <- ggplot() +
  geom_polygon(data = map.world, aes(x = long, y = lat, group = group), alpha = 0.5,  fill = c("grey")) +
  coord_fixed() +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, col = Order), shape = 21, size=1.5, col = c("#104E8B"), fill = c("#104E8B")) +
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
mapdf2 <- cbind(mapdf2[1], stack(mapdf2[2:6]))
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
  theme(axis.text.x = element_text(angle = 50, size = 8, vjust = 0.5, hjust = 0.7),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank())


# put map and bar chart together
scatter_bar_mercator <- cowplot::plot_grid(p_map1.1,p_bar, ncol = 1, labels = c("(a)", "(b)"))

# save plot (jpeg)
ggsave("FIGURE_1_scatter_bar_mercator_UN.jpeg", device ="jpeg", path = plotDir, width=175, height=225, units="mm", dpi = 350)



