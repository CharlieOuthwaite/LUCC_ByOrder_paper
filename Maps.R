##%######################################################%##
#                                                          #
####            Maps of site distribution               ####
#                                                          #
##%######################################################%##

# in this script, maps showing global site distribution and the proportional
# representation of each order in UN subregions are created

# clear working environment
rm(list = ls())

# set directories 
inDir<- "1_CheckPrepareData/"
outDir <- "1_CheckPrepareData/"
if(!dir.exists(outDir)) dir.create(outDir)

# load packages
packages_plot <- c("patchwork","dplyr","ggplot2","scatterpie","gt","broom.mixed","MASS","webshot","sjPlot")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# load data set
sites <- readRDS(file = paste0(outDir,"PREDICTSSiteData.rds"))

## plot map of site distribution: scatterplot ##

# plot the raster in ggplot

map.world <- map_data('world')

# map of sites, coloured by order
p_map <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  coord_fixed() +
  labs(x = "Longitude", y = "Latitude", fill = "Order") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, color = Order), shape = 20, size=1) +
  scale_color_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 7),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size =7))

#save map
ggsave("map_sites.jpeg", device ="jpeg", path = outDir, width=25, height=15, units="cm", dpi = 350)

## plot map of site by proportion: scatterpie ##

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
            long = mean(Longitude), lat = mean(Latitude), # calculate the mean coordinates, for plotting the piecharts
            radius = 25*sqrt(length(UN_subregion)/nrow(sites))) # calculate the radius for each subregion

# need to correct the position of the pie charts in europe to prevent overlap
# southern europe - move north and east
mapdf$lat[mapdf$UN_subregion == "Southern Europe"] <- 43
mapdf$long[mapdf$UN_subregion == "Southern Europe"] <- -5
# northern europe - move north and east
mapdf$lat[mapdf$UN_subregion == "Northern Europe"] <- 68
mapdf$long[mapdf$UN_subregion == "Northern Europe"] <- 4

# map showing proportion of each order by UN subregion
p_map2 <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  geom_scatterpie(aes(x=long, y=lat, group=UN_subregion, r=radius, fill=Order),
                  data=mapdf, cols=c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"), color=NA, alpha=.8) +
  scale_fill_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
  coord_fixed() +
  labs(x = "Longitude", y = "Latitude", fill = "Order") +
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 7),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size =7))

# save map
ggsave("map_proportionate.jpeg", device ="jpeg", path = outDir, width=25, height=15, units="cm", dpi = 350)

## layer maps of sites and by proportion: scatterplot + scatterpie ##
p_map3 <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  coord_fixed() +
  labs(x = "Longitude", y = "Latitude", fill = "Order") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude), shape = 20, size=1) +
  geom_scatterpie(aes(x=long, y=lat, group=UN_subregion, r=radius),
                  data=mapdf, cols=c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"), color=NA, alpha=.8) +
  scale_fill_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

# save maps (jpeg)
ggsave("maps_layered.jpeg", device ="jpeg", path = outDir, width=25, height=10, units="cm", dpi = 350)
