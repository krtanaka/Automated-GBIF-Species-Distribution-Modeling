require(tidyverse)
require(rgbif)
require(sf)
require(raster)
require(geodata)
require(ENMeval)
require(ecospat)
require(rJava)
require(dismo)
require(countrycode)
library(readr)
library(colorRamps)
library(ggmap)
library(ggspatial)
library(doParallel)

# rm(list = ls())

source("script/GBIF_SDM_Functions.R") # REPLACE W/ YOUR PATH TO GBIF_SDM_Functions.R

species_list <- c(
  "Unomia stolonifera",
  "Lutjanus gibbus",
  "Heniochus diphreutes",
  "Herklotsichthys quadrimaculatus",
  "Acropora globiceps",
  "Isopora crateriformis"
)

occ_df = read_csv("data/occurances_multi.csv") %>% 
  filter(Scientific.Name %in% species_list[1]) %>%
  as.data.frame()

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

load("/mnt/ldrive/ktanaka/env_rs.RData")
env_rs_i = env_rs
env_rs_i[["Bathymetry.Min"]][ env_rs_i[["Bathymetry.Min"]] <= -1000] <- NA
env_rs_i <- crop(env_rs_i, extent(floor(range(occ_df$Longitude)), floor(range(occ_df$Latitude))))
env_rs_i = mask(rast(env_rs_i), rast(env_rs_i[["Bathymetry.Min"]]))
env_rs_i = stack(env_rs_i)
plot(env_rs_i[[2]], col = matlab.like(100))

# ---- 6: Batch Run MaxEnt Models on all species ----
maxent_results = run_maxent(occ_df, env_rs_i)
save(maxent_results, file = "maxent_results.rda")

# Batch clip prediction extents
clipped_rasters_list <- spp_clip_raster(occ_df, env_rs, "Oahu", 500); plot(clipped_rasters_list[[1]])
# env_rs <- terra::resample(rast(env_rs), rast(bathy))

plot(maxent_results$models$`Lutjanus gibbus`)
plot(maxent_results$models$`Unomia stolonifera`)

r <- predict(maxent_results$models$`Lutjanus gibbus`, clipped_rasters_list[[1]])
r <- predict(maxent_results$models$`Unomia stolonifera`, clipped_rasters_list[[1]]) 

plot(r, col = matlab.like(100))
r = rasterToPoints(raster(r)) %>% as.data.frame()

# use ggmap
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Get the coordinates of the cell centers
coords <- coordinates(clipped_rasters_list[[1]] %>% stack())

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
                     shape = 22, alpha = 0.7, crs = 4326) + 
  scale_fill_gradientn(colors = matlab.like(100), "Habitat \nSuitability \n(0-1)") + 
  scale_color_gradientn(colors = matlab.like(100), "Habitat \nSuitability \n(0-1)") + 
  # ggtitle("Spatial distribution of U. stolonifera predicted habitat suitability") + 
  theme(legend.position = c(0.92, 0.81),
        legend.background = element_blank(), # Makes the legend background transparent
        legend.box.background = element_blank(), # Makes the legend box background transparent
        legend.text = element_text(color = "white"), # Makes the legend text white
        legend.title = element_text(color = "white") # Makes the legend title white
  )

ggsave(last_plot(), filename =  file.path("output/SDM_output.png"), height = 5.5, width = 5.5)

# Loop through each model and predict
maxent_predict()
