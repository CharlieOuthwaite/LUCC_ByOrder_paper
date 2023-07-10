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
dataDir<- "Data/"
inDir<- "1_CheckPrepareData/"
plotDir <- "1_CheckPrepareData/Plots/"
if(!dir.exists(plotDir)) dir.create(plotDir)


# load packages
packages_plot <- c("patchwork","dplyr","ggplot2","scatterpie","gt","broom.mixed","MASS","webshot","sjPlot","stringr","rjson","rgdal","maptools","sf","rgdal")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# load data set
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds"))

# load map
map.world <- map_data('world')

## Still trying to work out how to merge the polygons for subregions ###

# Checking the SpatialPolygon
# load map of country borders using rgdal
map.borders1 <- readOGR(dsn = paste0(dataDir,"/World Map/"), #Provide the directory of shp file
                        layer = "TM_WORLD_BORDERS-0.3", #Provide the name of the shp file without extension (.shp)
                        verbose = FALSE)

# fortify for use in ggplot (might be unnecessary step)
library(broom)
map.borders_df <- fortify(map.borders1)

# not working

# # code from first (abandoned) method, but keeping this here for now, just in case
# # try to merge country polygons based on UN subregion ID
# mb.coords <- coordinates(map.borders1)
# mb.id <- cut(mb.coords[,1], quantile(mb.coords[,1]), include.lowest=TRUE)
# map.borders1 <- maptools::unionSpatialPolygons(map.borders1, mb.id)

# second method

# load data as dataframe with package sf
map.borders <- sf::st_read(dsn = paste0(dataDir,"/World Map/"), layer = "TM_WORLD_BORDERS-0.3")

# try to merge country multipolygons based on value SUBREGION

# using st_union
# map.merge<-map.borders %>% 
#   group_by(SUBREGION) %>%
#   summarise(geometry = sf::st_union(geometry)) %>%
#   ungroup()

# using st_combine
map.merge <- map.borders %>% 
  group_by(SUBREGION) %>%
  summarise(geometry = sf::st_combine(geometry)) %>%
  ungroup()

# Simple feature collection with 23 features and 1 field
# Geometry type: MULTIPOLYGON
# Dimension:     XY
# Bounding box:  xmin: -180 ymin: -90 xmax: 180 ymax: 83.6236
# Geodetic CRS:  WGS 84
# # A tibble: 23 × 2
# SUBREGION                                                                                geometry
# <int>                                                                      <MULTIPOLYGON [°]>
#   1         0 (((96.91415 -12.19805, 96.90248 -12.2, 96.9147 -12.15195, 96.92165 -12.16806, 96.924...
#   2         5 (((-68.60861 -54.8914, -68.62056 -54.8914, -68.64311 -54.88861, -68.63722 -54.885, -...
#   3        11 (((2.484418 6.340486, 2.48 6.338611, 2.455 6.333055, 2.404722 6.33, 2.369722 6.33027...
#   4        13 (((-87.80334 17.29417, -87.80917 17.28917, -87.82195 17.29111, -87.82501 17.29166, -...
#   5        14 (((29.2299 -3.750964, 29.2325 -3.743333, 29.23638 -3.719722, 29.23889 -3.701945, 29....
#   6        15 (((2.96361 36.80222, 2.981389 36.80694, 3.001111 36.80971, 3.024167 36.80805, 3.0313...
#   7        17 (((11.75083 -16.75528, 11.775 -16.80473, 11.77 -16.8025, 11.75889 -16.79556, 11.7458...
#   8        18 (((37.85083 -46.95695, 37.84638 -46.96084, 37.83555 -46.96722, 37.82916 -46.96945, 3...
#   9        21 (((-64.85584 32.27861, -64.86111 32.27638, -64.87306 32.27694, -64.87556 32.28139, -...
#   10        29 (((-61.68667 17.02444, -61.73806 16.98972, -61.82917 16.99694, -61.87611 17.01694, -...
#   # … with 13 more rows

# add field with subregion names (found here: https://unstats.un.org/unsd/methodology/m49/)
SUBREGION <- c()

# add fields with lat-long coords
LAT <- c()
LONG <- c()


# plot to test
ggplot() +
  geom_polygon(data = map.borders_df, aes(x = long, y = lat,), fill="#69b3a2", color="white") +
  theme_void() 

# not working as of yet

## then incorporate into other maps as "geom_polygon()"


## 1. plot map of site distribution: scatterplot ##

# plot the raster in ggplot

# map of sites, no colour
p_map1.1 <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  coord_fixed() +
  labs(fill = "Order") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude), shape = 20, size=1) +
  #scale_color_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

#save map (mercator)
ggsave("map_sites_mercator.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)

# # change map projection to mollweide
# if (require("maps")) {
#   p_map1.1 + coord_map("mollweide",xlim=c(-180,180))
# }
# 
# #save map (mollweide)
# ggsave("map_sites_mollweide.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)
# 
# 
# # map of sites, coloured by order
# p_map1.2 <-ggplot() +
#   geom_map(data=map.world, map=map.world,
#            aes(x=long, y=lat, group=group, map_id=region),
#            fill= "grey", colour="grey", size=0.2) +
#   coord_fixed() +
#   labs(x = "Longitude", y = "Latitude", fill = "Order") +
#   geom_point(data = sites, aes(x = Longitude, y = Latitude, color = Order), shape = 20, size=1) +
#   scale_color_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
#   scale_x_continuous(breaks = c(-180,0,180), limits = c(-180,180)) +
#   scale_y_continuous(breaks = c(-90,-23.5,0,23.5,90), limits = c(-90,90)) +
#   theme(axis.title = element_text(size = 8),
#         axis.text.y = element_text(size =7),
#         axis.text.x = element_text(size =7),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         legend.position = "right",
#         legend.title = element_text(size = 8),
#         legend.text = element_text(size =7))
# 
# #save map (mercator)
# ggsave("map_sites_colour_mercator.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)
# 
# # change map projection to mollweide
# if (require("maps")) {
#   p_map1.2 + coord_map("mollweide",xlim=c(-180,180))
# }
# 
# #save map (mollweide)
# ggsave("map_sites_colour_mollweide.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)


## 2. plot map of sites by proportion: scatterpie ##
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
            long = mean(Longitude), lat = mean(Latitude), # calculate the mean coordinates, for plotting the piecharts
            radius = 25*sqrt(length(UN_subregion)/nrow(sites))) # calculate the radius for each subregion

# # need to correct the position of the pie charts in europe to prevent overlap
# # southern europe - move north and east
# mapdf$lat[mapdf$UN_subregion == "Southern Europe"] <- 43
# mapdf$long[mapdf$UN_subregion == "Southern Europe"] <- -5
# # northern europe - move north and east
# mapdf$lat[mapdf$UN_subregion == "Northern Europe"] <- 68
# mapdf$long[mapdf$UN_subregion == "Northern Europe"] <- 4
# 
# # map showing proportion of each order by UN subregion
# p_map2 <-ggplot() +
#   geom_map(data=map.world, map=map.world,
#            aes(x=long, y=lat, group=group, map_id=region),
#            fill= "grey", colour="grey", size=0.2) +
#   geom_scatterpie(aes(x=long, y=lat, group=UN_subregion, r=radius, fill=Order),
#                   data=mapdf, cols=c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"), color=NA, alpha=.8) +
#   scale_fill_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
#   scale_x_continuous(breaks = c(-180,0,180), limits = c(-180,180)) +
#   scale_y_continuous(breaks = c(-90,-23.5,0,23.5,90), limits = c(-90,90)) +
#   coord_fixed() +
#   labs(x = "Longitude", y = "Latitude", fill = "Order") +
#   theme(axis.title = element_text(size = 8), 
#         axis.text = element_text(size = 7),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         legend.position = "bottom",
#         legend.key.size = unit(3,"mm"),
#         legend.title = element_blank(),
#         legend.text = element_text(size =7))
# 
# # save map
# ggsave("map_proportionate_mercator.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)
# 
# # change map projection to mollweide
# if (require("maps")) {
#   p_map2 + coord_map("mollweide",xlim=c(-180,180))
# }
# 
# # save map
# ggsave("map_proportionate_mollweide.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)
# 
# 
# ## 3. layer maps of sites and by proportion: scatterplot + scatterpie ##
# p_map3 <-ggplot() +
#   geom_map(data=map.world, map=map.world,
#           aes(x=long, y=lat, group=group, map_id=region),
#            fill= "grey", colour="grey", size=0.2) +
#   coord_fixed() +
#   labs(x = "Longitude", y = "Latitude", fill = "Order") +
#   geom_point(data = sites, aes(x = Longitude, y = Latitude), shape = 20, size=1) +
#   geom_scatterpie(aes(x=long, y=lat, group=UN_subregion, r=radius),
#                   data=mapdf, cols=c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"), color=NA, alpha=.8) +
#   scale_fill_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5")) +
#   scale_x_continuous(breaks = c(-180,0,180), limits = c(-180,180)) +
#   scale_y_continuous(breaks = c(-90,-23.5,0,23.5,90), limits = c(-90,90)) +
#   theme_void() +
#   theme(axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(angle = 90, size = 8),
#         axis.text = element_blank(),
#         plot.background = element_blank(),
#         panel.background = element_blank(),
#         legend.position = 'bottom',
#         legend.box.spacing = unit(5,"mm"),
#         legend.key.size = unit(3,"mm"),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 7))
# 
# # save map
# ggsave("maps_layered_mercator.jpeg", device ="jpeg", path = plotDir, width=20, height=10, units="cm", dpi = 350)
# 
# if (require("maps")) {
#   p_map3 + coord_map("mollweide",xlim=c(-180,180))
# }
# 
# # save map
# ggsave("maps_layered_mollweide.jpeg", device ="jpeg", path = plotDir, width=20, height=10, units="cm", dpi = 350)

## 4. Scatterplot with barplot ##

# need to reorganize dataframe for barplot
mapdf2 <- mapdf %>% dplyr::select(UN_subregion, Coleoptera, Diptera, Hemiptera, Hymenoptera, Lepidoptera)
mapdf2 <- cbind(mapdf2[1], stack(mapdf2[2:3:4:5:6]))
mapdf2 <- mapdf2 %>% rename_at('ind', ~'Order')
mapdf2 <- dplyr::select(mapdf2,UN_subregion, Order, values)

mapdf2$UN_subregion <- factor(mapdf2$UN_subregion)
mapdf2$Order <- factor(mapdf2$Order) 

p_bar <- ggplot(data=mapdf2, aes(x=UN_subregion, y=values, fill=Order)) +
  scale_fill_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
  geom_bar(stat="identity")+
  labs(x = "UN subregion", y = "Number of sites") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(0, 1600), breaks = c(0, 250, 500, 750, 1000, 1250, 1500))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, size = 8),
        legend.position = "bottom",
        #legend.key.size = unit(3,"mm"),
        legend.title = element_blank()
        #legend.text = element_text(size = 8)
        )

# put map and bar chart together
scatter_bar_mercator <- cowplot::plot_grid(p_map1.1,p_bar, ncol = 1)

# save plot (jpeg)
ggsave("scatter_bar_mercator.jpeg", device ="jpeg", path = plotDir, width=300, height=300, units="mm", dpi = 350)

# # change projection
# p_map1.1 <- if (require("maps")) {
#   p_map1.1+ coord_map("mollweide",xlim=c(-180,180))
# }
# 
# # put map and bar chart together
# scatter_bar_mollweide <- cowplot::plot_grid(p_map1.1,p_bar, ncol = 1)
# 
# # save plot (jpeg)
# ggsave("scatter_bar_mollweide.jpeg", device ="jpeg", path = plotDir, width=300, height=300, units="mm", dpi = 350)
