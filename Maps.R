## Maps ##

packages_plot <- c("patchwork","dplyr","ggplot2","scatterpie","gt","broom.mixed","MASS","webshot","sjPlot",)
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# plot map of site distribution: scatterplot

# plot the raster in ggplot
map.world <- map_data('world')

# map of sites, coloured by order
p_map <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  coord_fixed() +
  labs(x = "Longitude", 
       y = "Latitude",
       fill = "Order") +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, colour = factor(Order)), shape = 20) +
  theme(axis.title = element_text(), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(),
      # legend.title = element_text()
        legend.position = "none") +
  ggtitle("a")

ggsave("map_sites.jpeg", device ="jpeg", path = outDir, width=25, height=15, units="cm", dpi = 350)

# plot map of site by proportion: scatterpie

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
            radius = 75*length(UN_subregion)/nrow(sites)) # calculate the radius for each subregion

# need to correct the position of the pie charts in europe
# southern europe - move north and east
mapdf$lat[mapdf$UN_subregion == "Southern Europe"] <- 43
mapdf$long[mapdf$UN_subregion == "Southern Europe"] <- -5
# northern europe - move north and east
mapdf$lat[mapdf$UN_subregion == "Northern Europe"] <- 68
mapdf$long[mapdf$UN_subregion == "Northern Europe"] <- 4

p_map2 <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  geom_scatterpie(aes(x=long, y=lat, group=UN_subregion, r=radius),
                  data=mapdf, cols=c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera"), color=NA, alpha=.8) +
  coord_fixed() +
  labs(x = "Longitude", 
       y = "Latitude",
       fill = "Order") +
  theme(axis.title = element_text(), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(),
      # legend.title = element_text()
        legend.position = "none") +
  ggtitle("b")

ggsave("map_proportionate.jpeg", device ="jpeg", path = outDir, width=25, height=15, units="cm", dpi = 350)

# get the legend
legend <- get_legend(
  p_map +
    guides(color = guide_legend(title = "Order", override.aes = list(size = 2, shape = 16)))+
    theme(legend.position = "left",
          legend.background = element_blank(), 
          legend.text = element_text(size = 10), 
          legend.title = element_text())
  )

# plot together
plot_figure <- plot_grid(
  plot_grid(p_map, p_map2, ncol = 1, rel_heights = c(1,1)),
  plot_grid(NULL,legend,NULL, ncol = 1, rel_heights = c(0.5,1,0.5)),
  nrow=1, rel_widths = c(1,0.2))

ggsave(filename = paste0(outDir, "Figure1_map_simplemods.pdf"), plot = last_plot(), width = 250, height = 200, units = "mm", dpi = 300)

# save plots
ggsave("maps.jpeg", device ="jpeg", path = outDir, width=25, height=20, units="cm", dpi = 350)