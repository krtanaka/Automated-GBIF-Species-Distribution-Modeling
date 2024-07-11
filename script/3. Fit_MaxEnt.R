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

# Check how many occurrences subset for each spp.
table(occ_df$Scientific.Name)

load("/mnt/ldrive/ktanaka/env_rs_i.RData")
load("/mnt/ldrive/ktanaka/env_rs.RData")

env_rs_i = env_rs
env_rs_i[["Bathymetry.Min"]][ env_rs_i[["Bathymetry.Min"]] <= -1000] <- NA
env_rs_i <- crop(env_rs_i, extent(range(occ_df$Longitude) + c(-1, 1), range(occ_df$Latitude) + c(-1, 1)))
env_rs_i = mask(rast(env_rs_i), rast(env_rs_i[["Bathymetry.Min"]]))
env_rs_i = stack(env_rs_i)
plot(env_rs_i[[3]], col = matlab.like(100))

# Batch Run MaxEnt Models on all species
maxent_results = run_maxent(occ_df, env_rs_i)
