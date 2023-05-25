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

# Load required libraries
packages <- c("predictsFunctions","patchwork", "dplyr", "ggplot2", "yarg", "lme4", "gt", "broom.mixed", "MASS","webshot")
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
# 934845 obs. of 67 variables

# Merge sites that have the same coordinates, from the same study and same taxonomic family (e.g. multiple traps on a single transect)
predicts <- predictsFunctions::MergeSites(diversity = predicts)
# 826292 obs. of 67 variables

# remove entries without Order
predicts <- droplevels(predicts[(predicts$Order!=""),])
# 826016 obs. of 67 variables

# check
table(predicts$Order)

# recode Isoptera to Blattodea and Phthiraptera to Psocodea
# Isoptera = infraorder within Blattodea; Phtiraptera = parvorder within Psocodea (Catalogue of Life)

predicts$Order <- recode_factor(predicts$Order,
                                Isoptera = "Blattodea",
                                Phthiraptera = "Psocodea")
# 826016 obs. of 67 variables

# convert Order to a "factor"
predicts$Order <- as.factor(predicts$Order)

#### 2. Calculate site metrics and prepare dataset ####

# Split predicts into separate data frames according to insect Order 

# use split function to split the predicts data frame into 6 data frames (1/Order)
OrderName <- paste0("",predicts$Order)

by_Order <- split(predicts,OrderName)

# extract data frames from list into global environment
list2env(by_Order,globalenv())

# Calculate site metrics of diversity for each order, include extra columns:
Archaeognatha <- droplevels(Archaeognatha)
Archaeognatha <- SiteMetrics(diversity = Archaeognatha,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Blattodea <- droplevels(Blattodea)
Blattodea <- SiteMetrics(diversity = Blattodea,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Coleoptera <- droplevels(Coleoptera)
Coleoptera <- SiteMetrics(diversity = Coleoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Dermaptera <- droplevels(Dermaptera)
Dermaptera <- SiteMetrics(diversity = Dermaptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Diptera <- droplevels(Diptera)
Diptera <- SiteMetrics(diversity = Diptera,
                       extra.cols = c("Predominant_land_use",
                                      "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Embioptera <- droplevels(Embioptera)
Embioptera <- SiteMetrics(diversity = Embioptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Ephemeroptera <- droplevels(Ephemeroptera)
Ephemeroptera <- SiteMetrics(diversity = Ephemeroptera,
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
Mantodea <- droplevels(Mantodea)
Mantodea <- SiteMetrics(diversity = Mantodea,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Mecoptera <- droplevels(Mecoptera)
Mecoptera <- SiteMetrics(diversity = Mecoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Neuroptera <- droplevels(Neuroptera)
Neuroptera <- SiteMetrics(diversity = Neuroptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Odonata <- droplevels(Odonata)
Odonata <- SiteMetrics(diversity = Odonata,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Orthoptera <- droplevels(Orthoptera)
Orthoptera <- SiteMetrics(diversity = Orthoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Phasmida <- droplevels(Phasmida)
Phasmida <- SiteMetrics(diversity = Phasmida,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Psocodea <- droplevels(Psocodea)
Psocodea <- SiteMetrics(diversity = Psocodea,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Siphonaptera <- droplevels(Siphonaptera)
Siphonaptera <- SiteMetrics(diversity = Siphonaptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Thysanoptera <- droplevels(Thysanoptera)
Thysanoptera <- SiteMetrics(diversity = Thysanoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Trichoptera <- droplevels(Trichoptera)
Trichoptera <- SiteMetrics(diversity = Trichoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Zoraptera <- droplevels(Zoraptera)
Zoraptera <- SiteMetrics(diversity = Zoraptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))
Zygentoma <- droplevels(Zygentoma)
Zygentoma <- SiteMetrics(diversity = Zygentoma,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order","UN_subregion","Best_guess_binomial"))

# merge all sites_Order data frames into one called "sites"
# merge using rbind()
sites <- rbind(Archaeognatha,Blattodea,Coleoptera,Dermaptera,Diptera,Embioptera,
               Ephemeroptera,Hemiptera,Hymenoptera,Lepidoptera,Mantodea,
               Mecoptera,Neuroptera,Odonata,Orthoptera,Phasmida,Psocodea,
               Siphonaptera,Thysanoptera,Trichoptera,Zoraptera,Zygentoma)


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

# 15346 obs. of 29 variables


# ...not sure, but it looks like we are re-classifying all secondary vegetation as "Light use"
# but why after we have already re-coded the land use type and intensity at all sites?
# sites$Use_intensity[((sites$LandUse=="Mature secondary vegetation") & 
#                        (sites$Use_intensity=="Intense use"))] <- "Light use"
# sites$Use_intensity[((sites$LandUse=="Intermediate secondary vegetation") & 
#                        (sites$Use_intensity=="Intense use"))] <- "Light use"
# sites$Use_intensity[((sites$LandUse=="Young secondary vegetation") & 
#                        (sites$Use_intensity=="Intense use"))] <- "Light use"

# remove the urban sites and sites that are NA in LUI
sites <- sites[!sites$LUI == "Urban", ]
sites <- sites[!is.na(sites$LUI), ]

# 12968 obs. of 29 variables

sites <- droplevels(sites)

# transform abundance values 
sites$LogAbund <- log(sites$Total_abundance+1)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]

# 12962 obs. of 30 variables

# create a new variable designating sites as Tropical or Non-tropical
# assign new variable for tropical/temperate, convert to factor, and filter out NA
sites$Realm <- ifelse(sites$Latitude >= -23.5 & sites$Latitude <= 23.5, "Tropical", "NonTropical")
sites$Realm <- factor(sites$Realm, levels = c("NonTropical", "Tropical"))
sites <- sites %>%
  filter(!is.na(Realm))

#### 3a. Summarise data by order ####

sites_summary <- sites %>%
  group_by(Order) %>% 
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
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
            Abundance = sum(LogAbund>0) ,
            SpeciesRichness = sum(Species_richness>0)) %>%
  ungroup() %>%  
  arrange(desc(Sites))%>%
  gt() %>%
  tab_header(title = "Data spread across land-uses"
  ) %>%
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
gtsave(sites_summary,"sites_summary_allorders.png",path = outDir)

# save as csv
write.csv(sites_summary,file = paste0("sites_summary_allorders.csv"))

#### 3b. Summarize data by order and latitudinal realm ####

lat_summary <- sites %>%
  group_by(Order,Realm) %>%
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
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
            Abundance = sum(LogAbund>0) ,
            SpeciesRichness = sum(Species_richness>0)) %>%
  ungroup() %>%
  gt() %>%
  tab_header(title = "Data spread across land-uses, by latitudinal realm"
  ) %>%
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
    columns = c(Order,Realm,Sites,Studies, Species,sites_1,studies_1,sites_2,studies_2,sites_low,studies_low,sites_high,studies_high,Abundance,SpeciesRichness)
  )%>%
  cols_label(
    Order = "Order",
    Realm = "Realm",
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
gtsave(lat_summary,"lat_summary_allorders.png",path = outDir)

# save as csv
write.csv(lat_summary,file = paste0("lat_summary_allorders.csv"))

# keep top five orders (according to number of sites sampled)
sites <- sites %>% filter(Order %in% c("Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera", "Hemiptera")) %>% droplevels()

# 8884 obs. of 31 variables

# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds")) 

#### 3c. Summarise remaining data by order ####

sites_summary <- sites %>%
  group_by(Order) %>% 
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
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
            Abundance = sum(LogAbund>0) ,
            SpeciesRichness = sum(Species_richness>0)) %>%
  ungroup() %>%  
  arrange(desc(Sites))%>%
  gt() %>%
  tab_header(title = "Data spread across land-uses"
  ) %>%
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
gtsave(sites_summary,"sites_summary_BigFive.png",path = outDir)

# save as csv
write.csv(sites_summary,file = paste0(outDir,"sites_summary_BigFive.csv"))

#### 3d. Summarize remaining data by order and latitudinal realm ####

lat_sites_summary <- sites %>%
  group_by(Order,Realm) %>%
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
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
            Abundance = sum(LogAbund>0) ,
            SpeciesRichness = sum(Species_richness>0)) %>%
  ungroup() %>%
  gt() %>%
  tab_header(title = "Data spread across land-uses, by latitudinal realm"
  ) %>%
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
    columns = c(Order,Realm,Sites,Studies,Species,sites_1,studies_1,sites_2,studies_2,sites_low,studies_low,sites_high,studies_high,Abundance,SpeciesRichness)
  )%>%
  cols_label(
    Order = "Order",
    Realm = "Realm",
    Sites = "Sites",
    Species = "Species",
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
gtsave(lat_sites_summary,"lat_sites_summary_BigFive.png", path = outDir)

# save as csv
write.csv(lat_sites_summary, file = paste0(outDir, "lat_sites_summary_BigFive.csv"))


#### site level. big five summaries ####

length(unique(sites$SS))
# 254 studies

length(unique(sites$SSBS))
# 6014 sites

pred_sub <- predicts[predicts$SSBS %in% sites$SSBS & predicts$Order %in% sites$Order, ]
pred_sub <- droplevels(pred_sub)
nrow(pred_sub) # 741933 records
table(pred_sub$Order)
# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#     377237       24149       45217      144474      150856
table(pred_sub$Predominant_land_use)
length(unique(pred_sub$SS)) # 254
length(unique(pred_sub$SSBS)) # 6014
length(unique(pred_sub$Species)) # 4324


length(unique(pred_sub$Taxon_name_entered)) # 17218

# determine unique taxa counts
species <- unique(pred_sub[,c('Order',"Family", 'Taxon_name_entered')])

# order level counts
order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#       6037        1477        1302        4528        3874


# family level counts
family.counts <- tapply(X = species$Family,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

# Coleoptera     Diptera   Hemiptera Hymenoptera Lepidoptera 
#        103          76          68          54          60
