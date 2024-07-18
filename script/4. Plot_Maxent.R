library(readr)
library(raster)
library(terra)
library(dplyr)
library(colorRamps)
library(ggmap)
library(ggspatial)

# run "2.Prep_Prediction_Layers.R" first

source("script/GBIF_SDM_Functions.R")

species_list <- c(
  "Unomia stolonifera",
  "Lutjanus gibbus",
  "Heniochus diphreutes",
  "Herklotsichthys quadrimaculatus",
  "Acropora globiceps",
  "Isopora crateriformis"
)

occ_df = read_csv("data/occurances_multi.csv") %>% 
  filter(Scientific.Name %in% species_list) %>%
  as.data.frame()

load("output/maxent_result_Herklotsichthys quadrimaculatus.rda")
load("output/maxent_result_Heniochus diphreutes.rda")
load("output/maxent_result_Isopora crateriformis.rda")
load("output/maxent_result_Unomia stolonifera.rda")

plot(maxent_result$model)

clipped_raster <- spp_clip_raster(occ_df, env_rs, "Tutuila", 500)#; plot(clipped_raster[[1]])
# env_rs <- terra::resample(rast(env_rs), rast(bathy))

r <- predict(maxent_result$model, clipped_raster)
plot(r, col = matlab.like(100))
r = rasterToPoints(raster(r)) %>% as.data.frame()

# use ggmap
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Get the coordinates of the cell centers
coords <- coordinates(clipped_raster %>% stack())

# Calculate the mean latitude and longitude
mean_lat <- mean(coords[, 2], na.rm = TRUE)
mean_lon <- mean(coords[, 1], na.rm = TRUE)

map = ggmap::get_map(location = c(mean_lon, mean_lat),
                     maptype = "satellite",
                     # zoom = 7,
                     force = T)
ggmap(map) +
  geom_spatial_point(data = r, aes(x, y, fill = maxent, color = maxent), 
                     size = 10,
                     shape = 22, alpha = 0.8, crs = 4326) + 
  scale_fill_gradientn(colors = matlab.like(100), "Predicted \nOccupancy \n(0-1)") + 
  scale_color_gradientn(colors = matlab.like(100), "Predicted \nOccupancy \n(0-1)") + 
  # ggtitle("Spatial distribution of U. stolonifera predicted habitat suitability") + 
  theme(legend.position = c(0.92, 0.81),
        legend.background = element_blank(), # Makes the legend background transparent
        legend.box.background = element_blank(), # Makes the legend box background transparent
        legend.text = element_text(color = "white"), # Makes the legend text white
        legend.title = element_text(color = "white") # Makes the legend title white
  )

ggsave(last_plot(), filename =  file.path("output/SDM_output.png"), height = 5.5, width = 5.5)

