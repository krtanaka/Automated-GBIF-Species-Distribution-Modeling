# https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html
# https://data-blog.gbif.org/post/gbif-filtering-guide/

# install.packages("usethis")
# usethis::edit_r_environ()

library(rgbif)
library(dplyr)
library(CoordinateCleaner)
library(ggplot2)

taxonkey <- name_backbone("Unomia stolonifera")$usageKey
taxonkey <- name_backbone("Lutjanus gibbus")$usageKey
taxonkey <- name_backbone("Heniochus diphreutes")$usageKey
taxonkey <- name_backbone("Herklotsichthys quadrimaculatus")$usageKey
taxonkey <- name_backbone("Acropora globiceps")$usageKey
taxonkey <- name_backbone("Isopora crateriformis")$usageKey

# set up gbif credentials first
# https://docs.ropensci.org/rgbif/articles/gbif_credentials.html

gbif_download <- occ_download(
  pred("taxonKey", taxonkey),
  pred("hasCoordinate", TRUE), 
  pred("hasGeospatialIssue", FALSE), # remove GBIF default geospatial issues
  format = "SIMPLE_CSV") 

occ_download_wait(gbif_download) 

# filtering pipeline  
df <- gbif_download %>%
  occ_download_get() %>%
  occ_download_import() %>%
  setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
  filter(occurrencestatus  == "PRESENT") %>%
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")) %>%
  # filter(year >= 2000) %>%
  filter(year >= 1900) %>%
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
  filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>%
  # cc_cen(buffer = 2000) %>% # remove country centroids within 2km
  # cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
  # cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
  # cc_sea() %>% # remove from ocean 
  distinct(decimallongitude, decimallatitude, specieskey, datasetkey, countrycode, .keep_all = TRUE) %>%
  glimpse() # look at results of pipeline

species = unique(df$species)
species = gsub(" ", "_", species)

df <- df %>% 
  dplyr::rename(Scientific.Name = species,
                Longitude  = decimallongitude, 
                Latitude  = decimallatitude,
                Country = countrycode) %>%
  dplyr::select(Scientific.Name, Longitude, Latitude, Country) %>%
  na.omit()

df %>% 
  mutate(Longitude = ifelse(Longitude < 0, Longitude + 360, Longitude)) %>% 
  ggplot(aes(Longitude, Latitude)) + 
  geom_point(size = 5, color = "blue", alpha = 0.8) + 
  annotation_map(map_data("world")) + 
  coord_quickmap()

readr::write_csv(df, file = paste0("occurances_", species, ".csv"))
