# ---- Run MaxEnt Batch ----

run_maxent = function(occ_sf, env) {
  
  set.seed(2024)
  
  # occ_sf = occ_df
  # env = env_rs
  # env = env_rs_i
  
  # Get unique Scientific Names
  print(unique(occ_sf$Scientific.Name))
  species = unique(occ_sf$Scientific.Name)
  
  # # Initialize lists to store results
  # enm_results = list()
  # maxent_models = list()
  
  # Loop through each species
  for(sp in species) {
    
    # sp = "Unomia stolonifera"
    # sp = species[1]
    
    cat(paste("Processing", sp, "\n"))
    
    # Subset occurrence data for the species
    occ_sp = occ_sf[occ_sf$Scientific.Name == sp, c("Longitude", "Latitude")]
    
    # Run ENMevaluate
    enmeval_results = ENMevaluate(occ_sp, env, 
                                  bg = NULL, 
                                  tune.args = list(fc = c("L","LQ","H", "LQH", "LQHP", "LQHPT"), rm = 1:5), 
                                  partitions = "randomkfold", 
                                  partition.settings = list(kfolds = 2), 
                                  algorithm = "maxnet", 
                                  n.bg = 100,
                                  parallel = T,
                                  numCores = detectCores()/2,
                                  updateProgress = T,
                                  taxon.name = sp)  # specify the taxon name here
    
    enmeval_df = enmeval_results@results
    
    # Store ENMevaluate results
    # enm_results[[sp]] = enmeval_df
    
    # Subset the ENMeval results to get the best model
    enmeval_bestm = subset(enmeval_df, delta.AICc == 0)
    
    # Decode the features
    maxent_feats = as.character(enmeval_bestm$fc)
    maxent_rm = as.character(enmeval_bestm$rm)
    
    # Print out the results
    cat(paste("Best features for", sp, ":", maxent_feats, "\n"))
    cat(paste("Best regularization multiplier for", sp, ":", maxent_rm, "\n"))
    
    # Run MaXent SDM
    sp_maxent_model = dismo::maxent(env, as.matrix(occ_sp), features = maxent_feats, betamultiplier = maxent_rm)
    
    # Store maxent model
    # maxent_models[[sp]] = sp_maxent_model
    maxent_result = list(enm = enmeval_df, model = sp_maxent_model)
    save(maxent_result, file = paste0("output/maxent_result_", sp, ".rda"))
    
  }
  
  # return(list(enm = enm_results, models = maxent_models))
}

# ---- Clip Prediction Raster by Species ----

spp_clip_raster <- function(spp, env_rs, domain, depth_cutoff) {
  
  domain <- "Tutuila"
  
  box <- readr::read_csv("data/Bounding_Boxes.csv")
  
  box_i <- box %>% 
    group_by(unit) %>% 
    summarise(ymin = min(ymin), 
              ymax = max(ymax),
              xmin = min(xmin),
              xmax = max(xmax)) 
  
  box_r <- box %>% 
    group_by(region) %>% 
    summarise(ymin = min(ymin), 
              ymax = max(ymax),
              xmin = min(xmin),
              xmax = max(xmax)) %>% 
    rename(unit = region)
  
  box <- rbind(box_i, box_r) %>% filter(unit == domain)
  
  # Define the extent from island_boxes
  extent_box <- extent(box$xmin, box$xmax, box$ymin, box$ymax)
  
  # Clip env_rs using the defined extent
  clipped_raster <- crop(env_rs, extent_box)
  
  clipped_raster[["Bathymetry.Min"]][ clipped_raster[["Bathymetry.Min"]] <= -depth_cutoff] <- NA
  clipped_raster = mask(rast(clipped_raster), rast(clipped_raster[["Bathymetry.Min"]]))
  
  return(clipped_raster)
}
