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

setwd("~/Automated-GBIF-Species-Distribution-Modeling/")
setwd("C:/Users/kisei.tanaka/Automated-GBIF-Species-Distribution-Modeling")

# ---- 1: Load Script Functions ----

source("GBIF_SDM_Functions.R")        # REPLACE W/ YOUR PATH TO GBIF_SDM_Functions.R

# ---- 2: Download BioClim Data ----

bioclim_download(resolution = 2.5) # Adjust resolution as needed

# ---- 3: Batch Read BioClim Folder Data ----

bioclim_dir = here::here("data2", "climate", "wc2.1_2.5m") # Adjust version & resolution path accordingly
bioclim_dir = here::here("data2", "climate", "wc2.1_5m") # Adjust version & resolution path accordingly
bioclim_dir = here::here("data2", "climate", "bio_oracle_v3") # Adjust version & resolution path accordingly

bioclim_files = list.files(path = bioclim_dir, pattern = "\\.tif$", full.names = TRUE)
bioclim_files = list.files(path = bioclim_dir, pattern = "\\.nc$", full.names = TRUE)

# Convert to RasterStack for dismo::maxent()
env_rs = raster::stack(bioclim_files)
env_rs <- crop(env_rs, extent(-160, 125, -5, 24))
env_rs[[5]][env_rs[[5]] <= -100] <- NA
env_rs[[1]] <- mask(env_rs[[1]], env_rs[[5]])
env_rs[[2]] <- mask(env_rs[[2]], env_rs[[5]])
env_rs[[3]] <- mask(env_rs[[3]], env_rs[[5]])
env_rs[[4]] <- mask(env_rs[[4]], env_rs[[5]])
env_rs[[5]] <- mask(env_rs[[5]], env_rs[[5]])
env_rs[[6]] <- mask(env_rs[[6]], env_rs[[5]])

# bathy = raster("~/data2/ETOPO_2022_v1_15s_3f38_8816_d630.nc")
# bathy[bathy <= -30] <- NA
# bathy[bathy >= 0] <- NA  
# bathy = readAll(bathy)
# save(bathy, file = "~/Automated-GBIF-Species-Distribution-Modeling/data2/etopo_0-30.rdata")
load("~/Automated-GBIF-Species-Distribution-Modeling/data2/etopo_0-30.rdata")
env_rs <- crop(env_rs, extent(bathy))
env_rs <- resample(env_rs, bathy)
env_rs <- mask(env_rs, bathy)
save(env_rs, file = "~/Automated-GBIF-Species-Distribution-Modeling/data2/env_rs.rdata")
load("~/Automated-GBIF-Species-Distribution-Modeling/data2/env_rs.rdata")

# Optionally, save the clipped raster
writeRaster(clipped_env_rs, "clipped_env_rs.grd", overwrite = TRUE)

# ---- 4: Load Species Dataframe & Convert Countries ----

# Example Dataframe - Replace W/ Your Own
spp = data.frame(
  # Scientific.Name = "Adelges tsugae",
  Scientific.Name = "Unomia stolonifera",
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

# ---- 6: Batch Run MaxEnt Models on all species ----

maxent_results = run_maxent(occ_df, env_rs)
save(maxent_results, file = "maxent_results.rda")

# ---- 7: Predict Species in Source.Location ----

# Use world country boundaries to clip prediction extent for each species
# https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/information/?flg=en-us
worldbound = st_read(here::here('data2', 'world-administrative-boundaries', 'world-administrative-boundaries.shp'))

# Batch clip prediction extents
clipped_rasters_list <- spp_clip_raster(spp, worldbound, env_rs)
clipped_rasters_list <- spp_clip_raster_island(spp, worldbound, env_rs, "Oahu")

# Create the "MaxEnt_Predictions" directory to store results
dir.create("MaxEnt_Target_Predictions", showWarnings = FALSE)

# Loop through each model and predict
maxent_predict()

