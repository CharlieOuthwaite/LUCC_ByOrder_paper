##%######################################################%##
#                                                          #
#### Exploration of data from primary vegetation sites  ####
#                                                          #
##%######################################################%##


# Charlie Outhwaite 25/10/22

# In this script, I will look into the spread of climate data available for 
# sites in primary vegetation. This is to check whether some studies with more
# broadly distributed sites are driving the negative/positive trends in this 
# land use. 

# ensure environment is clear
rm(list = ls())

# load libraries
library(ggplot2)

# directories
datadir <- "5_RunLUIClimateModels/"
outdir <- "11_Explore_primveg/"
if(!dir.exists(outdir)) dir.create(outdir)

# load in the dataset that includes the climate info
pred <- readRDS(file = paste0(datadir,"PREDICTSSitesClimate_Data.rds"))
pred <- droplevels(pred)

# number of unique studies
length(unique(pred$SS)) # 257

# number of sites within studies
site_tab <- data.frame(table(pred$SS))
range(site_tab$Freq) # number of sites ranges from 1 to 680


#### plot sites within a study ####


# for each study, plot sites

# plot the raster in ggplot
map.world <- map_data('world')

# i <- site_tab$Var1[1]

for(i in site_tab$Var1){
  
  # subset the data
  datsub <- pred[pred$SS == i & pred$Predominant_land_use == "Primary vegetation", ]
  
  if(nrow(datsub) == 0) next
  
  
  # create and save a map
  p1 <-ggplot() +
    geom_map(data=map.world, map=map.world,
             aes(x=long, y=lat, group=group, map_id=region),
             fill= "grey", colour="grey", size=0.2) +
    geom_point(data = datsub, aes(x = Longitude, y = Latitude, colour = factor(Order)), shape = 20, size = 2) +
    theme_bw() + 
    theme(axis.title = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank(), 
          legend.title = element_blank())
  
  ggsave(p1, filename = paste0(outdir, i, "pointmap.png"))
  
}



#### look at range of anomaly values for each study. ####



# subset to studies with sites in primary vegetation
prim <- pred[pred$Predominant_land_use == "Primary vegetation", ]
prim <- droplevels(prim)

library(dplyr)

prim_sum <- prim %>% group_by(SS) %>% summarise(min = min(StdTmeanAnomaly), max = max(StdTmeanAnomaly), dif = (max(StdTmeanAnomaly- min(StdTmeanAnomaly))) )

View(prim_sum)

# there are two studies which have a very large difference across site anomaly values
# AD1_2010__Davis 1 (sites spread across the UK), SC1_2006__Benedick 1 (some sites in Borneo)

View(prim[prim$SS %in% c("AD1_2010__Davis 1", "SC1_2006__Benedick 1"), ]) # 22 sites


# suggest a sensitivity test where these two are removed. 
