library(ggplot2)
library(ggmap)
library(maps)
library(ggthemes)
library(tidyverse)
library(ggrepel)

setwd("~/Desktop/dros_speciation/manuscript/figures/making_figures/figure1")

data <- read.csv("Matute_subset_populations.csv")

#Subset data for labels
data$Label <- ifelse(data$Sample %in% c("SD", "SP", "ZH", "ZS", "ZI", "LZV", "MC", "DRM", "ZW"), "", data$Sample)

world <- map_data("world")
world <- subset(world, lat > -60)

countries <- world %>%
  group_by(region) %>%
  summarize(
    Longitude = mean(range(long)),
    Latitude = mean(range(lat))
  )



bbox <- data.frame(
  xmin = 20,  # Minimum longitude
  xmax = 40,  # Maximum longitude
  ymin = -30, # Minimum latitude
  ymax = -13  # Maximum latitude
)

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "antiquewhite", color = "lightgray") +
  geom_point(data = data, aes(x = Longitude, y = Latitude), shape=20, color="orange3", size = 3) +
  geom_rect(data = bbox, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = NA, color = "black", size = 0.5) +
  geom_text(data = data, aes(x = Longitude, y = Latitude, label = Label), vjust = 0, hjust = 1.1, size = 5, color = "black",
            nudge_y = ifelse(data$Label == "EF", 1.5, 0),
            nudge_x = ifelse(data$Label == "EF", 5, 0)) +
  theme_map() + theme(plot.background = element_rect(fill = "lightblue3", color = NA),
                      legend.position = "none")






ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "antiquewhite", color = "lightgray") +
  geom_point(data = data, aes(x = Longitude, y = Latitude), shape=20, color="orange3", size = 3) +
  geom_text(data = data, aes(x = Longitude, y = Latitude, label = Sample), vjust = 0, hjust = 1.1, size = 5, color = "black",
            nudge_y = ifelse(data$Sample == "DRM", 0.1, ifelse(data$Sample == "MC", -0.4, 0)),
            nudge_x = ifelse(data$Sample == "DRM", 1.6, ifelse(data$Sample == "MC", 1.2, 0))) +
  coord_quickmap(
    xlim = c(20, 40),  # Longitude limits
    ylim = c(-30, -13) # Latitude limits
  ) +
  geom_text(data = countries, aes(x = Longitude, y = Latitude, label = region), size = 3, color = "green4",
            nudge_x = ifelse(countries$region == "Malawi", -0.3, ifelse(countries$region == "Mozambique", -0.3, ifelse(countries$region == "South Africa", -3, 0))),
            nudge_y = ifelse(countries$region == "Zambia", -0.7, ifelse(countries$region == "South Africa", 6, 0))) +
  theme_map() + theme(plot.background = element_rect(fill = "lightblue3", color = NA),   # Plot background
                      legend.position = "none")





