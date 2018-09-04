## extract temperature variables for sites with butterfly data

## housekeeping
library(raster); library(sp); library(lubridate); library(dplyr); library(reshape2)
library(ggplot2)
UK_outline <- readRDS("files/UK_outline.rds")

# get coordinates to extract temperature data at
df <- readRDS("files/clean butterfly dataset_v2.rds")
siteKey <- df %>%
    dplyr::select(siteID, lon_BNG, lat_BNG) %>%
    unique
CRS_BNG <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 
               +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
site_BNG <- siteKey %>%
    dplyr::select(-siteID) %>%
    SpatialPoints(., proj4string = CRS_BNG)


## (1) extract climatic data ----
# (1.1) get file names ----
temp_names <- paste0("../UKCP09/temp_", seq(1970, 2010, 10), "s/")
fnames_temp <- lapply(temp_names, list.files, full.names=T) %>% unlist
years <- gsub(".*temperature_(....).*", "\\1", fnames_temp)

# (1.2) extract temp data ----
extracted_temp <- list()
for(i in seq(fnames_temp)) {
    print(years[i])
    r <- brick(fnames_temp[i])
    proj4string(r) <- CRS_BNG
    extracted <- raster::extract(r, site_BNG)
    extracted_temp[[i]] <- as_data_frame(extracted) %>%
        mutate(siteID = 1:n()) %>%
        reshape2::melt(., id.vars="siteID") %>%
        as_data_frame %>%
        group_by(siteID) %>%
        mutate(day = 1:length(siteID), year=years[i])
}

# format
UKCP_temp <- bind_rows(extracted_temp) %>%
    ungroup %>% 
    mutate(siteID = siteKey$siteID[siteID],
           date = as.Date(paste(year, day), format="%Y %j"),
           month = month(date)) %>%
    dplyr::select(siteID, year, month, day, date, temp=value)

## (2) create summary df ----
temp_summ <- UKCP_temp %>% 
    group_by(siteID, year, month) %>% summarise(mean=mean(temp), sd=sd(temp), min = min(temp), max=max(temp))

## (3) wide-format df----
sd_wide <- temp_summ %>% 
    # remove mean (i.e. work with sd)
    dplyr::select(-mean) %>%
    dcast(siteID + year ~ month) %>%
    rename_at(as.character(1:12), function(x) paste0("sd_", x)) %>%
    mutate(year=as.numeric(year))

extra_sd <- sd_wide %>%
    dplyr::select(siteID:sd_8) %>%
    mutate(year=as.numeric(year),
           year = year - 1)
names(extra_sd)[3:10] <- paste0("sd_", 13:20)

sd_verywide <- sd_wide %>%
    left_join(., extra_sd) %>% 
    as_data_frame

# now mean
mean_wide <- temp_summ %>% 
    # remove mean (i.e. work with sd)
    dplyr::select(-sd) %>%
    dcast(siteID + year ~ month) %>%
    rename_at(as.character(1:12), function(x) paste0("temp_", x)) %>%
    mutate(year=as.numeric(year))

extra_mean <- mean_wide %>%
    dplyr::select(siteID:temp_8) %>%
    mutate(year=as.numeric(year),
           year = year - 1)
names(extra_mean)[3:10] <- paste0("temp_", 13:20)

mean_verywide <- mean_wide %>%
    left_join(., extra_mean) %>% 
    as_data_frame

temp_full <- left_join(mean_verywide, sd_verywide) %>% 
    as_data_frame

## full monthly df
saveRDS(temp_full, "files/temp_wideformat.rds")
# save
saveRDS(temp_summ, "files/monthly temperature.rds")

## (4) sense checks ----
# does extracted temperature look sensible? 
# spatial & temporal plots
temp_site <- inner_join(siteKey, temp_summ) %>%
    mutate(year = as.integer(year))

temp_site %>%
    mutate(period = cut(year, seq(1970, 2020, 10), include.lowest=T, dig.lab = 5)) %>%
    group_by(siteID, lon_BNG, lat_BNG, period) %>%
    summarise(meanT = mean(mean)) %>%
    ggplot(aes(lon_BNG, lat_BNG, col=meanT)) + 
    geom_polygon(data = UK_outline, aes(x = lon_BNG, y = lat_BNG, group = group), fill="grey80", col="black") +
    geom_point() +
    facet_wrap(~period) +
    coord_fixed() +
    scale_color_gradient2(midpoint=10, low = "darkblue", high="indianred") +
    labs(caption=bquote(bold(Figure)~"Average decadal temperatures for retained sites"), 
         x="Easting", y="Northing", colour=bquote("Mean"~degree*C)) +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          strip.text=element_text(hjust=0, face="bold"))

ggsave("figures/temperature checks/spatial plot of decadal temperatures.png", 
       width=210, height=180, units="mm")

set.seed(101)
temp_site %>%
    filter(siteID %in% sample(unique(siteID), 10), year %in% sample(unique(year), 20)) %>%
    mutate(siteID = paste0("siteID: ", siteID)) %>%
    # group_by(siteID, lon_BNG, lat_BNG, period) %>%
    # summarise(meanT = mean(mean)) %>%
    ggplot(aes(month, mean, col=factor(year))) + 
    geom_line() + 
    facet_wrap(~siteID) +
    theme_minimal() +
    guides(col=F) +
    theme(strip.text=element_text(hjust=0, face="bold")) +
    scale_x_continuous(breaks=seq(1, 12, 1)) +
    labs(caption=bquote(bold(Figure)~"Monthly mean temperatures for 10 sites and 20 years in dataset"), 
         x="Month", y=bquote("Monthly mean"~degree*C)) 
ggsave("figures/temperature checks/temporal plot of monthly temperatures.png", 
       width=210, height=180, units="mm")   
