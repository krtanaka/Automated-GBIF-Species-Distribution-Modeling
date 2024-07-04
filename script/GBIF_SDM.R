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

rm(list = ls())

source("GBIF_SDM_Functions.R") # REPLACE W/ YOUR PATH TO GBIF_SDM_Functions.R

occ_df = read_csv("occurances_Unomia_stolonifera.csv") %>% as.data.frame()
occ_df = read_csv("occurances_Heniochus_diphreutes.csv") %>% as.data.frame()

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

# ---- 3: Batch Read Folder Data ----
bio_oracle_dir = file.path(fs::path_home(), "Desktop/data", "bio_oracle_v3") 
bio_oracle_files = list.files(path = bio_oracle_dir, pattern = "\\.nc$", full.names = TRUE)
bio_oracle_rs = raster::stack(bio_oracle_files)

anth_dir = file.path(fs::path_home(), "Desktop/data", "sedac_gfw") 
anth_files = list.files(path = anth_dir, full.names = TRUE)
anth_rs <- list()
pb <- txtProgressBar(min = 0, max = length(anth_files), style = 3)
for (e in 1:length(anth_files[1:2])) {
  
  # e = 2
  
  setTxtProgressBar(pb, e)
  anth_rs_e <- terra::rast(anth_files[e])
  
  if (e %in% c(3:6)) {
    
    name = names(anth_rs_e[[1]])
    anth_rs_e = mean(anth_rs_e[[1:5]])
    names(anth_rs_e) = name
    
  }
  
  anth_rs_e <- terra::resample(anth_rs_e, rast(bio_oracle_rs))
  anth_rs[[e]] <- raster(anth_rs_e)
  
}
close(pb)
anth_rs <- stack(anth_rs)

env_rs = stack(bio_oracle_rs, anth_rs)

env_rs_i = env_rs
env_rs_i[["Bathymetry.Min"]][ env_rs_i[["Bathymetry.Min"]] <= -100] <- NA
env_rs_i <- crop(env_rs_i, extent(floor(range(occ_df$Longitude)), floor(range(occ_df$Latitude))))

pb <- txtProgressBar(min = 0, max = dim(env_rs_i)[3], style = 3)

for (e in 1:dim(env_rs_i)[3]) {
  
  # e = 2
  
  setTxtProgressBar(pb, e)
  env_rs_i[[e]] <- mask(rast(env_rs_i[[e]]), rast(env_rs_i[["Bathymetry.Min"]]))
  # plot(env_rs_i[[e]], col = matlab.like(100))
  
}

close(pb)
plot(env_rs_i, col = matlab.like(100))

# ---- 6: Batch Run MaxEnt Models on all species ----
maxent_results = run_maxent(occ_df, env_rs_i)
save(maxent_results, file = "maxent_results.rda")


# Use world country boundaries to clip prediction extent for each species
# https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/information/?flg=en-us
worldbound = st_read(file.path(fs::path_home(), "Desktop/data", 'world-administrative-boundaries', 'world-administrative-boundaries.shp'))

# Batch clip prediction extents
# clipped_rasters_list <- spp_clip_raster(spp, worldbound, env_rs)
clipped_rasters_list <- spp_clip_raster_island(occ_df, worldbound, env_rs, "Oahu")
plot(clipped_rasters_list[[1]])
# env_rs <- terra::resample(rast(env_rs), rast(bathy))

# Create the "MaxEnt_Predictions" directory to store results
dir.create("MaxEnt_Target_Predictions", showWarnings = FALSE)

plot(maxent_results$models[1]$`Unomia stolonifera`)

r <- predict(maxent_results$models[1]$`Unomia stolonifera`, clipped_rasters_list[[1]]) 
plot(r, col = matlab.like(100))

# use ggmap
ggmap::register_google("AIzaSyDpirvA5gB7bmbEbwB1Pk__6jiV4SXAEcY")

# Get the coordinates of the cell centers
coords <- coordinates(clipped_rasters_list[[1]])

# Calculate the mean latitude and longitude
mean_lat <- mean(coords[, 2], na.rm = TRUE)
mean_lon <- mean(coords[, 1], na.rm = TRUE)

map = ggmap::get_map(location = c(mean_lon, mean_lat),
                     maptype = "satellite",
                     zoom = 9,
                     force = T)

r = rasterToPoints(r) %>% as.data.frame()

ggmap(map) +
  geom_spatial_point(data = r, aes(x, y, fill = layer, color = layer), 
                     size = 8, shape = 22, alpha = 0.7, crs = 4326) + 
  scale_fill_gradientn(colors = matlab.like(100), "Habitat \nSuitability \n(0-1)") + 
  scale_color_gradientn(colors = matlab.like(100), "Habitat \nSuitability \n(0-1)") + 
  ggtitle("Spatial distribution of U. stolonifera predicted habitat suitability") + 
  theme(legend.position = c(0.92, 0.81),
        legend.background = element_blank(), # Makes the legend background transparent
        legend.box.background = element_blank(), # Makes the legend box background transparent
        legend.text = element_text(color = "white"), # Makes the legend text white
        legend.title = element_text(color = "white") # Makes the legend title white
  )

ggsave(last_plot(), filename =  file.path(fs::path_home(), "Desktop/Unomia_Thermal_SDM_output.png"), height = 8, width = 8)

# Loop through each model and predict
maxent_predict()



