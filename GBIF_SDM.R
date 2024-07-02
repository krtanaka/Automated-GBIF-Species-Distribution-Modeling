# ---- 0: Load Packages ----

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

rm(list = ls())

# ---- 1: Load Script Functions ----

source("GBIF_SDM_Functions.R") # REPLACE W/ YOUR PATH TO GBIF_SDM_Functions.R

# ---- 4: Load Species Dataframe & Convert Countries ----

# Example Dataframe - Replace W/ Your Own
spp = data.frame(
  Scientific.Name = "Unomia stolonifera",
  # Scientific.Name = "Herklotsichthys quadrimaculatus",
  Source.Location = c("Venezuela", "Cuba", "Philippines", "Indonesia", "Chinese Taipei")
)

# Convert Country --> Country Code
spp = spp %>%
  filter(Source.Location != "")

spp$countryCode = sapply(spp$Source.Location, flexible_country_code)

spp = subset(spp, !is.na(countryCode) & countryCode != "")

# Apply the function to the Source.Location column
spp$Source.Location = sapply(spp$Source.Location, region_to_country, mapping = region_mapping)

# ---- 5: Download Species Data from GBIF ----

# Apply the function to each scientific name in the dataframe
occ_list = setNames(mapply(gbif_occ_data, 
                           spp$Scientific.Name, 
                           spp$countryCode, 
                           SIMPLIFY = FALSE), 
                    spp$Scientific.Name)

# Convert the list to a data frame, and also create a new column to hold the Scientific.Name
occ_df = bind_rows(occ_list, .id = "Scientific.Name")

# Remove NA values
occ_df = occ_df[complete.cases(occ_df), ]

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

# ---- 3: Batch Read BioClim Folder Data ----
# Adjust version & resolution path accordingly
bioclim_dir = file.path(fs::path_home(), "Desktop/data", "bio_oracle_v3") 
sedac_gfw_dir = file.path(fs::path_home(), "Desktop/data", "sedac_gfw") 

bioclim_files = list.files(path = bioclim_dir, pattern = "\\.nc$", full.names = TRUE)
sedac_gfw_files = list.files(path = sedac_gfw_dir, full.names = TRUE)

# Convert to RasterStack for dismo::maxent()
env_rs = raster::stack(bioclim_files)
env_rs <- crop(env_rs, extent(floor(range(occ_df$Longitude)), floor(range(occ_df$Latitude))))
env_rs[["Bathymetry.Min"]][ env_rs[["Bathymetry.Min"]] <= -100] <- NA
for (e in 1:length(bioclim_files)) {
  
  env_rs[[e]] <- mask(env_rs[[e]], env_rs[["Bathymetry.Min"]])

}

anth_rs <- list()

for (e in 1:length(sedac_gfw_files)) {
  
  # e = 4
  
  anth_rs_e = terra::rast(sedac_gfw_files[e])
  anth_rs_e <- crop(anth_rs_e, extent(floor(range(occ_df$Longitude)), floor(range(occ_df$Latitude))))
  anth_rs_e = terra::resample(anth_rs_e, rast(env_rs))
  anth_rs_e <- mask(anth_rs_e, rast(env_rs[["Bathymetry.Min"]]))
  plot(anth_rs_e)
  
  anth_rs[[e]] <- raster(anth_rs_e)

}

anth_rs <- stack(anth_rs)
plot(anth_rs, col = matlab.like(100))

env_rs = stack(env_rs, anth_rs)

names(env_rs)

plot(env_rs[], col = matlab.like(100))

# bathy = raster("~/data2/ETOPO_2022_v1_15s_3f38_8816_d630.nc")
# bathy[bathy <= -30] <- NA
# bathy[bathy >= 0] <- NA  
# bathy = readAll(bathy)
# save(bathy, file = "~/Automated-GBIF-Species-Distribution-Modeling/data2/etopo_0-30.rdata")
# load(file.path(fs::path_home(), "Desktop/data/etopo_0-30.rdata"))
# bathy = raster(file.path(fs::path_home(), "Desktop/data/etopo360_3a6a_0a57_385e.nc"))

# env_rs <- crop(env_rs, extent(bathy))
# env_rs <- terra::resample(rast(env_rs), rast(bathy))
# env_rs = raster(env_rs)
# env_rs <- mask(env_rs, bathy)
# save(env_rs, file = "~/Automated-GBIF-Species-Distribution-Modeling/data2/env_rs.rdata")
# load("~/Automated-GBIF-Species-Distribution-Modeling/data2/env_rs.rdata")

# ---- 6: Batch Run MaxEnt Models on all species ----

maxent_results = run_maxent(occ_df, env_rs)
save(maxent_results, file = "maxent_results.rda")

# ---- 7: Predict Species in Source.Location ----

# Use world country boundaries to clip prediction extent for each species
# https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/information/?flg=en-us
worldbound = st_read(file.path(fs::path_home(), "Desktop/data", 'world-administrative-boundaries', 'world-administrative-boundaries.shp'))

# Convert to RasterStack for dismo::maxent()
env_rs = raster::stack(bioclim_files)
env_rs[["Bathymetry.Min"]][ env_rs[["Bathymetry.Min"]] <= -500] <- NA

for (e in 1:length(bioclim_files)) {
  
  env_rs[[e]] <- mask(env_rs[[e]], env_rs[["Bathymetry.Min"]])
  
}

# plot(env_rs)

# Batch clip prediction extents
# clipped_rasters_list <- spp_clip_raster(spp, worldbound, env_rs)
clipped_rasters_list <- spp_clip_raster_island(spp, worldbound, env_rs, "Hawaii")
plot(clipped_rasters_list[[1]])
# env_rs <- terra::resample(rast(env_rs), rast(bathy))

# Create the "MaxEnt_Predictions" directory to store results
dir.create("MaxEnt_Target_Predictions", showWarnings = FALSE)

plot(maxent_results$models[1]$`Unomia stolonifera`)

r <- predict(maxent_results$models[1]$`Unomia stolonifera`, clipped_rasters_list[[1]]) 
# plot(r, col = matlab.like(100))

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
                     size = 10, shape = 22, alpha = 0.8, crs = 4326) + 
  scale_fill_gradientn(colors = matlab.like(100), "Habitat \nSuitability \n(0-1)") + 
  scale_color_gradientn(colors = matlab.like(100), "Habitat \nSuitability \n(0-1)") + 
  ggtitle("Spatial distribution of U. stolonifera predicted habitat suitability") + 
  theme(legend.position = c(0.92, 0.81),
        legend.background = element_blank(), # Makes the legend background transparent
        legend.box.background = element_blank(), # Makes the legend box background transparent
        legend.text = element_text(color = "white"), # Makes the legend text white
        legend.title = element_text(color = "white") # Makes the legend title white
        )

ggsave(last_plot(), filename = '~/Desktop/Unomia_Thermal_SDM_output.png', height = 8, width = 8)

# Loop through each model and predict
maxent_predict()

