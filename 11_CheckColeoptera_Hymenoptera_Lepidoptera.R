##%######################################################%##
#                                                          #
####          Check Coleoptera Hymenoptera,                #
#                 and Lepidoptera Data                  ####
#                                                          #
##%######################################################%##

# This script checks the family composition of Hymenoptera and Lepidoptera

# ensure working directory is clear
rm(list = ls())

# set up directories
dataDir <- "C:/Users/Kyra/Documents/GitHub/LUCC_ByOrder_paper/Data/"
outDir <- "C:/Users/Kyra/Documents/GitHub/LUCC_ByOrder_paper/12_CheckColeoptera_Hymenoptera_Lepidoptera/"
if(!dir.exists(outDir)) dir.create(outDir)

# Load required libraries
packages <- c("predictsFunctions","patchwork", "dplyr", "tidyr", "ggplot", "yarg", "lme4", "gt", "broom.mixed", "MASS","webshot")
suppressWarnings(suppressMessages(lapply(packages, require, character.only = TRUE)))

# source in additional functions
source("C:/Users/Kyra/Documents/GitHub/LUCC_ByOrder_paper/0_Functions.R")

#### 1. Load and check dataset ####

# Set the path to copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)
# 3250404 obs. of 67 variables

# Select only data for insects
predicts <- predicts[(predicts$Class=="Insecta"),]
# 935078 obs. of 67 variables

predicts <- predicts %>% filter(Order %in% c("Coleoptera","Hymenoptera", "Lepidoptera")) %>% droplevels()

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
predicts <- MergeSites(diversity = predicts)
# 826292 obs. of 67 variables

# remove entries without Order
predicts <- droplevels(predicts[(predicts$Order!=""),])
# 826016 obs. of 67 variables

# convert Order to a "factor"
predicts$Order <- as.factor(predicts$Order)
predicts$Family <- as.factor(predicts$Family)

#### 2. Calculate site metrics ####

# Split predicts into separate data frames according to insect Order 

# use split function to split the predicts data frame into 2 data frames (1/Order)
OrderName <- paste0("",predicts$Order)

by_Order <- split(predicts,OrderName)

# extract data frames from list into global environment
list2env(by_Order,globalenv())

# Calculate site metrics of diversity for each order, include extra columns:
Coleoptera <- droplevels(Coleoptera)
Coleoptera <- SiteMetrics(diversity = Coleoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order","Family","Species"))
Hymenoptera <- droplevels(Hymenoptera)
Hymenoptera <- SiteMetrics(diversity = Hymenoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order","Family","Species"))
Lepidoptera <- droplevels(Lepidoptera)
Lepidoptera <- SiteMetrics(diversity = Lepidoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order","Family","Species"))

# merge all sites_Order data frames into one called "sites"
# merge using rbind()
sites <- rbind(Coleoptera,Hymenoptera,Lepidoptera)

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

# 15346 obs. of 26 variables


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

# 12968 obs. of 26 variables

sites <- droplevels(sites)

# transform abundance values 
sites$LogAbund <- log(sites$Total_abundance+1)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]

# 12962 obs. of 27 variables

# create a new variable designating sites as Tropical or Non-tropical
# assign new variable for tropical/temperate, convert to factor, and filter out NA
sites$Realm <- ifelse(sites$Latitude >= -23.5 & sites$Latitude <= 23.5, "Tropical", "NonTropical")
sites$Realm <- factor(sites$Realm, levels = c("NonTropical", "Tropical"))
sites <- sites %>%
  filter(!is.na(Realm))

# rename missing values in Family
sites <- sites %>% mutate_at(c('Family'), ~na_if(., ''))
sites$Family <- as.character(sites$Family)
sites$Family <- sites$Family %>% replace_na('Unspecified')
sites$Family <- as.factor(sites$Family)

### 3. Check Coleopteraa ###

# filter for just coleoptera
coleoptera <- sites[(sites$Order=="Coleoptera"),]

coleoptera_summary <- coleoptera %>%
  group_by(Family) %>% 
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
  summarise(sites = length(Family),
            studies = n_distinct(SS),
            percent = 100*(length(Family)/nrow(coleoptera)),
            species = n_distinct(Species),
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
  arrange(desc(sites))%>%
  gt() %>%
  tab_header(title = "Coleoptera summary"
  ) %>%
  tab_spanner(
    label = "Total",
    columns = c(sites,studies,species,percent)
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
    columns = c(Family,sites,studies,percent,species,sites_1,studies_1,sites_2,studies_2,sites_low,studies_low,sites_high,studies_high,Abundance,SpeciesRichness)
  ) %>%
  cols_label(
    Family = "Family",
    sites = "Sites",
    percent = "%",
    species = "Species",
    studies = "Studies",
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
gtsave(coleoptera_summary,"coleoptera_summary.png")

# save as csv
write.csv(coleoptera_summary,"coleoptera_summary.csv", row.names = TRUE)

### 4. Check Hymenoptera ###

# filter for just hymenopterans
hymenoptera <- sites[(sites$Order=="Hymenoptera"),]

hymenoptera_summary <- hymenoptera %>%
  group_by(Family) %>% 
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
  summarise(sites = length(Family),
            studies = n_distinct(SS),
            percent = 100*(length(Family)/nrow(hymenoptera)),
            species = n_distinct(Species),
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
  arrange(desc(sites))%>%
  gt() %>%
  tab_header(title = "Hymenoptera summary"
  ) %>%
  tab_spanner(
    label = "Total",
    columns = c(sites,studies,species,percent)
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
    columns = c(Family,sites,studies,percent,species,sites_1,studies_1,sites_2,studies_2,sites_low,studies_low,sites_high,studies_high,Abundance,SpeciesRichness)
  ) %>%
  cols_label(
    Family = "Family",
    sites = "Sites",
    percent = "%",
    species = "Species",
    studies = "Studies",
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
gtsave(hymenoptera_summary,"hymenoptera_summary.png")

# save as csv
write.csv(hymenoptera_summary,"hymenoptera_summary.csv", row.names = TRUE)

### 5. Check Lepidoptera ###

# filter for just lepidopterans
lepidoptera <- sites[(sites$Order=="Lepidoptera"),]

lepidoptera_summary <- lepidoptera %>%
  group_by(Family) %>% 
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
  summarise(sites = length(Family),
            studies = n_distinct(SS),
            percent = 100*(length(Family)/nrow(lepidoptera)),
            species = n_distinct(Species),
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
  arrange(desc(sites))%>%
  gt() %>%
  tab_header(title = "Lepidoptera summary"
  ) %>%
  tab_spanner(
    label = "Total",
    columns = c(sites,studies,species,percent)
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
    columns = c(Family,sites,studies,percent,species,sites_1,studies_1,sites_2,studies_2,sites_low,studies_low,sites_high,studies_high,Abundance,SpeciesRichness)
  ) %>%
  cols_label(
    Family = "Family",
    sites = "Sites",
    percent = "%",
    species = "Species",
    studies = "Studies",
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
gtsave(lepidoptera_summary,"lepidoptera_summary.png")

# save as csv
write.csv(lepidoptera_summary,"lepidoptera_summary.csv", row.names = TRUE)
