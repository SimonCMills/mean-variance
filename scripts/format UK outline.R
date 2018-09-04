# create dataframe of UK outline (for figures)

# packages
library(dplyr); library(ggplot2); library(sp)

# datum
crs_BNG <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 
               +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
crs_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# load mapdata and reformat to British National Grid
UK <- map_data(map = "world", region = "UK") %>%
    as_data_frame %>%
    rename(lon_WGS84 = long, lat_WGS84 = lat)
UK_BNG <- SpatialPoints(UK[,c(1,2)], crs_WGS84) %>%
    spTransform(., crs_BNG) %>%
    as_data_frame %>%
    select(lon_BNG = lon_WGS84, lat_BNG = lat_WGS84)

# bind together and save
UK_full <- bind_cols(UK_BNG, UK)
saveRDS(UK_full, "files/UK_outline.rds")
