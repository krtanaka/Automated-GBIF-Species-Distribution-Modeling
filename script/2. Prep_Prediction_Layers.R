library(terra)
library(raster)
library(usdm)

rm(list = ls())

bio_oracle_dir = file.path(fs::path_home(), "Desktop/data", "bio_oracle_v3") 
bio_oracle_files = list.files(path = bio_oracle_dir, pattern = "\\.nc$", full.names = TRUE)
bio_oracle_rs = raster::stack(bio_oracle_files); names(bio_oracle_rs)

v = vifstep(terra::rast(bio_oracle_rs), th = 10
            , keep = c(
              # "Maximum.OceanTemperature",
              # "Minimum.OceanTemperature",
              "Bathymetry.Min")
)

bio_oracle_rs = raster::subset(bio_oracle_rs, v@results$Variables); names(bio_oracle_rs)

anth_dir = file.path(fs::path_home(), "Desktop/data", "sedac_gfw") 
anth_files = list.files(path = anth_dir, full.names = TRUE)[1:2]
anth_rs <- list()

pb <- txtProgressBar(min = 0, max = length(anth_files), style = 3)
for (e in 1:length(anth_files)) {
  
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
# env_rs = readAll(env_rs)
# save(env_rs, file = "/Users/kisei.tanaka/Desktop/env_rs.RData")

