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
