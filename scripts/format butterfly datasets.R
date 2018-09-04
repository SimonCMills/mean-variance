## Format butterfly dataset
##

# housekeeping
library(xlsx); library(dplyr); library(ggplot2)
UK_outline <- readRDS("files/UK_outline.rds")

# read in raw data, site data key, and species data key
df <- read.csv("files/UKBMS/ukbms_sindex_2016_negative counts removed.csv") %>%
    as_data_frame %>%
    dplyr::select(species_code = SPECIES, 
                  brood = BROOD, 
                  siteID = SITE, 
                  year = YEAR, 
                  N_t = SINDEX)

siteKey <- read.table("files/UKBMS/UKBMS_site_list_updated_2013_2131sites.txt", header=T) %>%
    as_data_frame %>%
    dplyr::select(siteID = SITENO, 
                  lon_BNG = EAST, 
                  lat_BNG = NORTH, 
                  gridID = GRIDREF)

speciesKey <- read.xlsx2("files/UKBMS/species codes.xls", 1) %>%
    as_data_frame() %>%
    dplyr::select(species_code = species.code, 
                  species = common.name, 
                  species_LB = latin.name) %>%
    mutate(species_code = as.integer(as.character(species_code)), 
           species = as.character(species), 
           species_LB = as.character(species_LB))

# merge species key and site key into dataset
df_full <- df %>% 
    filter(brood == 0) %>%
    dplyr::select(-brood) %>%
    left_join(., speciesKey) %>%
    left_join(., siteKey)

## Quality check ----
# 6 species are dropped for not having a species name associated with the code
df_full %>% 
    filter(is.na(species)) %>%
    .$species_code %>% unique

# 851 sites are dropped for not having coordinates
df_full %>% 
    filter(is.na(lon_BNG)) %>%
    .$siteID %>% unique %>% length

# almost entirely sites with very few years of data; just 48 site:species combinations
# with more than 10 years of data (compared to 10270 initially)
summarised <- df_full %>% 
    filter(is.na(lon_BNG)) %>%
    group_by(siteID, species) %>%
    summarise(n_site = n())

summarised %>% ggplot(aes(n_site)) + geom_histogram()
summarised %>% filter(n_site >= 10) %>% nrow()

## Format dataframe ----
# function to get count in year prior (checks that year interval ==1)
count_prior <- function(count, diff_priorYr, interval=1) {
    out <- rep(NA, length(count))
    # out vector with an observation in year prior (when interval == 1) gets the 
    # count from the year prior. When there is no observation in the year prior, 
    # it returns NA
    out[which(diff_priorYr==interval)] <- count[(which(diff_priorYr==interval)-interval)]
    out
}

## format year prior in
df_fmtd <- df_full %>%
    filter(complete.cases(.)) %>%
    group_by(species, siteID) %>%
    arrange(species, siteID, year) %>%
    mutate(diff_priorYr = year - c(NA, year[-length(year)]),
           N_pre = count_prior(N_t, diff_priorYr, 1))

    
## apply cleaning step ----
# there have to be at least 10 functional observations (a non-zero count in either year
# t or year t-1) & not a long-distance migrant (3 cases)
df_cleaned <- df_fmtd %>%
    ungroup %>%
    na.omit %>% # remove NAs
    filter(N_t != 0 | N_pre != 0) %>% 
    group_by(siteID, species) %>%
    filter(n() >= 10) %>%
    group_by(species) %>%
    filter(length(unique(siteID)) >= 50) %>%
    ungroup %>%
    filter(!species %in% c("Clouded yellow", "Red admiral", "Painted lady"))

## plot layout & nDP summary ----
df_cleaned %>%
    dplyr::select(siteID, lon_BNG, lat_BNG) %>%
    unique %>%
    ggplot(aes(lon_BNG, lat_BNG)) +
    geom_polygon(data = UK_outline, aes(x = lon_BNG, y = lat_BNG, group = group), fill="white", col="black") +
    geom_point() +
    coord_equal() +
    labs(caption=bquote(bold(Figure)~"Sites retained after formatting, rough version"), 
         x="Easting", y="Northing") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank())
ggsave("figures/retained sites after formatting_rough ver.png")

df_cleaned %>%
    dplyr::select(species, siteID, lon_BNG, lat_BNG) %>%
    unique %>%
    mutate(id = as.numeric(as.factor(species))) %>%
    filter(id %in% 1:20) %>%
    ggplot(aes(lon_BNG, lat_BNG)) +
    geom_polygon(data = UK_outline, aes(x = lon_BNG, y = lat_BNG, group = group), fill="white", col="black") +
    geom_point() +
    coord_equal() +
    labs(caption=bquote(bold(Figure)~"Sites retained after formatting, rough version"), 
         x="Easting", y="Northing") +
    facet_wrap(~species, ncol=5) +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          strip.text=element_text(hjust=0, face="bold", size=8))

ggsave("figures/retained sites after formatting_by species_p1.png", 
       height=297, width=200, units="mm")

df_cleaned %>%
    dplyr::select(species, siteID, lon_BNG, lat_BNG) %>%
    unique %>%
    mutate(id = as.numeric(as.factor(species))) %>%
    filter(id %in% 21:32) %>%
    ggplot(aes(lon_BNG, lat_BNG)) +
    geom_polygon(data = UK_outline, aes(x = lon_BNG, y = lat_BNG, group = group), fill="white", col="black") +
    geom_point() +
    coord_equal() +
    labs(caption=bquote(bold(Figure)~"Sites retained after formatting, rough version"), 
         x="Easting", y="Northing") +
    facet_wrap(~species, ncol=5) +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          strip.text=element_text(hjust=0, face="bold", size=8))

ggsave("figures/retained sites after formatting_by species_p2.png", 
       height=297*.75, width=200, units="mm")


summarise_nDP <- df_cleaned %>%
    group_by(species) %>%
    summarise(`Number of observations` = n(), `Number of sites` = length(unique(siteID))) %>%
    arrange(`Number of observations`,`Number of sites`) %>% 
    filter(`Number of sites` > 50) %>%
    mutate(species=factor(species, levels=species))  %>%
    reshape2::melt()

ggplot(summarise_nDP, aes(value, species)) + geom_point() + facet_wrap(~variable, scale="free_x") +
    theme_bw() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(hjust=0, face="bold"))
ggsave("figures/nDP_summary.png", width=190, height=150, units="mm")

# create summary table
summaryTable <- df_cleaned %>%
    group_by(species) %>%
    summarise(fullName = unique(paste0(species, ", ", species_LB)), 
              `Number of observations` = n(), 
              `Number of sites` = length(unique(siteID)), 
              `first year` = min(year),
              `last year` = max(year)) %>%
    arrange(`Number of observations`,`Number of sites`) %>% 
    filter(`Number of sites` > 50) %>% 
    arrange(desc(`Number of observations`))
write.csv(summaryTable, "files/summary_table.csv")
## save----
saveRDS(df_cleaned, "files/clean butterfly dataset_v2.rds")
