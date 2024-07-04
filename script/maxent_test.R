
fnames = list.files(path = "/Users/kisei.tanaka/Desktop/BioOracle/", pattern = "\\.nc$", full.names = TRUE, recursive = T)
fnames <- fnames[grep("all_units", fnames)]

variable_names <- sapply(fnames, function(x) {
  basename <- gsub(".*/", "", x)
  gsub("_all_units\\.nc", "", basename)
})

# Print the result
variable_names

simplified_names <- gsub("_depthMean_Baseline_\\d{4}-\\d{4}_.*|_Baseline_\\d{4}-\\d{4}_.*", "", variable_names)
simplified_names <- gsub("SeaWaterDirection", "SWD", simplified_names)
simplified_names <- gsub("SeaWaterSpeed", "SWS", simplified_names)
simplified_names <- gsub("DissolvedMolecularOxygen", "Oxygen", simplified_names)
simplified_names

predictors <- stack(fnames)
names(predictors) = simplified_names
plot(predictors)

occ = read_csv("/Users/kisei.tanaka/Desktop/unomia_stolonifera_occurances.csv")
occ = occ %>% select(lon, lat) %>% as.data.frame()

me <- maxent(predictors, occ)
me
plot(me)

fnames = list.files(path = "/Users/kisei.tanaka/Desktop/BioOracle/", pattern = "\\.nc$", full.names = TRUE, recursive = T)
fnames <- fnames[grep("Unit_Level_Data/Oahu", fnames)]

plot(stack(stack(fnames[1]) %>% mean(),
     stack(fnames[2]) %>% mean(),
     stack(fnames[3]) %>% mean(),
     stack(fnames[4]) %>% mean(),
     stack(fnames[5]) %>% mean(),
     stack(fnames[6]) %>% mean(),
     stack(fnames[7]) %>% mean(),
     stack(fnames[8]) %>% mean(),
     stack(fnames[9]) %>% mean(),
     stack(fnames[10]) %>% mean(),
     stack(fnames[11]) %>% mean()))

predictors <- stack(stack(stack(fnames[1]) %>% mean(),
                              stack(fnames[2]) %>% mean(),
                              stack(fnames[3]) %>% mean(),
                              stack(fnames[4]) %>% mean(),
                              stack(fnames[5]) %>% mean(),
                              stack(fnames[6]) %>% mean(),
                              stack(fnames[7]) %>% mean(),
                              stack(fnames[8]) %>% mean(),
                              stack(fnames[9]) %>% mean(),
                              stack(fnames[10]) %>% mean(),
                              stack(fnames[11]) %>% mean()))
names(predictors) = simplified_names
plot(predictors)

r <- predict(me, predictors) 

# with some options:
# r <- predict(me, predictors, args=c("outputformat=raw"), progress='text', 
#      filename='maxent_prediction.grd')

plot(r, col = matlab.like(100))
points(occ)
