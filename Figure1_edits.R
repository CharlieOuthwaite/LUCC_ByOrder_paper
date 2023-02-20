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
plotDir <- "1_CheckPrepareData/Plots/"
if(!dir.exists(plotDir)) dir.create(plotDir)

# load packages
packages_plot <- c("patchwork","dplyr","ggplot2","scatterpie","gt","broom.mixed","MASS","webshot","sjPlot","stringr")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# load data set
sites <- readRDS(file = paste0(plotDir,"PREDICTSSiteData.rds"))

## 1. plot map of site distribution: scatterplot ##

# plot the raster in ggplot

map.world <- map_data('world')

# map of sites, no colour
p_map1.1 <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  coord_fixed() +
  labs(x = "Longitude", y = "Latitude", fill = "Order") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude), shape = 20, size=1) +
  #scale_color_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(angle = 90, size = 8),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

#save map (mercator)
ggsave("map_sites_mercator.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)

# change map projection to mollweide
if (require("maps")) {
  p_map1.1 + coord_map("mollweide",xlim=c(-180,180))
}

#save map (mollweide)
ggsave("map_sites_mollweide.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)


# map of sites, coloured by order
p_map1.2 <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  coord_fixed() +
  labs(x = "Longitude", y = "Latitude", fill = "Order") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, color = Order), shape = 20, size=1) +
  scale_color_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5"))+
  scale_x_continuous(breaks = c(-180,0,180), limits = c(-180,180)) +
  scale_y_continuous(breaks = c(-90,-23.5,0,23.5,90), limits = c(-90,90)) +
  theme(axis.title = element_text(size = 8),
        axis.text.y = element_text(size =7),
        axis.text.x = element_text(size =7),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size =7))

#save map (mercator)
ggsave("map_sites_colour_mercator.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)

# change map projection to mollweide
if (require("maps")) {
  p_map1.2 + coord_map("mollweide",xlim=c(-180,180))
}

#save map (mollweide)
ggsave("map_sites_colour_mollweide.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)


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
  scale_x_continuous(breaks = c(-180,0,180), limits = c(-180,180)) +
  scale_y_continuous(breaks = c(-90,-23.5,0,23.5,90), limits = c(-90,90)) +
  coord_fixed() +
  labs(x = "Longitude", y = "Latitude", fill = "Order") +
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 7),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(3,"mm"),
        legend.title = element_blank(),
        legend.text = element_text(size =7))

# save map
ggsave("map_proportionate_mercator.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)

# change map projection to mollweide
if (require("maps")) {
  p_map2 + coord_map("mollweide",xlim=c(-180,180))
}

# save map
ggsave("map_proportionate_mollweide.jpeg", device ="jpeg", path = plotDir, width=25, height=15, units="cm", dpi = 350)


## 3. layer maps of sites and by proportion: scatterplot + scatterpie ##
p_map3 <-ggplot() +
  geom_map(data=map.world, map=map.world,
          aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  coord_fixed() +
  labs(x = "Longitude", y = "Latitude", fill = "Order") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude), shape = 20, size=1) +
  geom_scatterpie(aes(x=long, y=lat, group=UN_subregion, r=radius),
                  data=mapdf, cols=c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"), color=NA, alpha=.8) +
  scale_fill_manual(values=c("#2271B2","#F0E442","#359B73","#E69F00","#F748A5")) +
  scale_x_continuous(breaks = c(-180,0,180), limits = c(-180,180)) +
  scale_y_continuous(breaks = c(-90,-23.5,0,23.5,90), limits = c(-90,90)) +
  theme_void() +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(angle = 90, size = 8),
        axis.text = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'bottom',
        legend.box.spacing = unit(5,"mm"),
        legend.key.size = unit(3,"mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7))

# save map
ggsave("maps_layered_mercator.jpeg", device ="jpeg", path = plotDir, width=20, height=10, units="cm", dpi = 350)

if (require("maps")) {
  p_map3 + coord_map("mollweide",xlim=c(-180,180))
}

# save map
ggsave("maps_layered_mollweide.jpeg", device ="jpeg", path = plotDir, width=20, height=10, units="cm", dpi = 350)

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
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(angle = 90, size = 8),
        axis.text.x = element_text(angle = 45, size = 8),
        # axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(3,"mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

# put map and bar chart together
scatter_bar_mercator <- cowplot::plot_grid(p_map1.1,p_bar, ncol = 1)

# save plot (jpeg)
ggsave("scatter_bar_mercator.jpeg", device ="jpeg", path = plotDir, width=300, height=300, units="mm", dpi = 350)

# change projection
p_map1.1 <- if (require("maps")) {
  p_map1.1+ coord_map("mollweide",xlim=c(-180,180))
}

# put map and bar chart together
scatter_bar_mollweide <- cowplot::plot_grid(p_map1.1,p_bar, ncol = 1)

# save plot (jpeg)
ggsave("scatter_bar_mollweide.jpeg", device ="jpeg", path = plotDir, width=300, height=300, units="mm", dpi = 350)
