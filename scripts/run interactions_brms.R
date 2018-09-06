# get model comparisons using brms and WAIC
#
#
library(brms); library(loo); library(dplyr)
options(mc.cores = parallel::detectCores())

dredged <- readRDS("../BFLY_spatial/files/dredgedCoefs_monthly_temp.rds") %>%
    filter(!grepl("Inter|log", coef))

# bfly_df <- readRDS("files/butterfly & temperature- following exclusion.rds")
bfly_df <- readRDS("files/clean butterfly dataset_v2.rds")
# check this is correct
bfly_df %>%
    dplyr::select(lon_BNG, lat_BNG) %>%
    unique %>%
    plot

# read monthly temperature data
month_df <- readRDS("files/temp_wideformat.rds")

alts <- readRDS("../BFLY_spatial/files/elevation data/elevations.rds")

bfly_full  <- month_df %>% 
    mutate(year = year-1) %>%
    left_join(bfly_df, ., by=c("siteID", "year")) %>% 
    inner_join(., alts)


## iterate
sp_list <- unique(bfly_full$species)
catch_waic <- list()
catch_preds <- list()
for(i in sp_list) {
    print(i)
    # get vars and extract data
    mean_i <- dredged %>% filter(species==i) %>% .$coef
    sd_i <- paste0(gsub("temp|rain", "sd", mean_i), "_", gsub("_[0-9]*", "", mean_i))
    CHB <- bfly_full %>%
        filter(species==i) %>%
        filter(elev_site <= 200) %>%
        dplyr::select(siteID, year, N_t, species, N_pre, lat_BNG, lon_BNG, elev_site,
                      mean = matches(paste0(mean_i, "$")), sd = matches(paste0(sd_i, "$"))) %>%
        na.omit %>% 
        mutate(siteID=factor(siteID), 
               mean = mean - mean(mean), 
               elev_site = scale(elev_site),
               lon_BNG = scale(lon_BNG), 
               lat_BNG = scale(lat_BNG),
               X_t =log(N_t+1), 
               X_pre = log(N_pre+1),
               X_pre_sc = scale(X_pre),
               X_diff = X_t - X_pre)
    
    b_lin_sd <- brm(X_diff~ X_pre_sc + mean*sd + 
                        # lat_BNG + mean:lat_BNG + lon_BNG + mean:lon_BNG + 
                        # elev_site + mean:elev_site + 
                        (1|siteID), iter=5000, data=CHB)
    b_lin <- brm(X_diff~ X_pre + mean + 
                     # lat_BNG + mean:lat_BNG + lon_BNG + mean:lon_BNG + 
                     # elev_site + mean:elev_site + 
                     (1|siteID), iter=5000, data=CHB)
    b_null <- brm(X_diff~ X_pre +
                      # lat_BNG + mean:lat_BNG + lon_BNG + mean:lon_BNG + 
                      # elev_site + mean:elev_site + 
                      (1|siteID), iter=5000, data=CHB)
    
    w_sd <- waic(b_lin_sd)
    w_noSD <- waic(b_lin)
    w_null <- waic(b_null)
    comp <- compare(w_sd, w_noSD)
    comp2 <- compare(w_sd, w_null)

    waic_i <- data_frame(species=i, dWAIC_1 = comp[1], dSE_1 = comp[2], dWAIC_2 = comp2[1], dSE_2 = comp2[2])
    predicts_i <- CHB %>% mutate(X_pre = median(X_pre))
    pSD_i <- fitted(b_lin_sd, re_form=NA, newdata=predicts_i) %>% as_data_frame
    pNOSD_i <- fitted(b_lin, re_form=NA, newdata=predicts_i) %>% as_data_frame
    names(pSD_i) <- paste0("sd_", c("est", "SE", "lwr", "upr"))
    names(pNOSD_i) <- paste0("mean_", c("est", "SE", "lwr", "upr"))
    predicts_full_i <- bind_cols(predicts_i, pSD_i, pNOSD_i)

    save(b_lin_sd, b_lin, b_null, w_sd, w_noSD, w_null, file = paste0("files/models/", i, "_controls.Rdata"))
    catch_waic[[i]] <- waic_i
    catch_preds[[i]] <- predicts_full_i
    
    saveRDS(catch_waic, "files/waic_updating.rds")
    saveRDS(catch_preds, "files/preds_updating.rds")
}


# Repeat, but with series of control variables
catch_waic <- list()
catch_preds <- list()
for(i in sp_list) {
    print(i)
    # get vars and extract data
    mean_i <- dredged %>% filter(species==i) %>% .$coef
    sd_i <- gsub("temp", "sd", mean_i) 
    CHB <- bfly_full %>%
        filter(species==i) %>%
        dplyr::select(siteID, year, N_t, species, N_pre, lat_BNG, lon_BNG, elev_site,
                      mean = matches(paste0(mean_i, "$")), sd = matches(paste0(sd_i, "$"))) %>%
        na.omit %>% 
        mutate(siteID=factor(siteID), 
               mean = mean - mean(mean), 
               elev_site = scale(elev_site),
               lon_BNG = scale(lon_BNG), 
               lat_BNG = scale(lat_BNG),
               X_t =log(N_t+1), 
               X_pre = log(N_pre+1),
               X_pre_sc = scale(X_pre),
               X_diff = X_t - X_pre)
    
    b_lin_sd <- brm(X_diff~ X_pre_sc + mean + sd + mean:sd + 
                        lat_BNG + mean:lat_BNG + lon_BNG + mean:lon_BNG +
                        elev_site + mean:elev_site +
                        (1|siteID), data=CHB)
    b_lin <- brm(X_diff~ X_pre + mean + 
                     lat_BNG + mean:lat_BNG + lon_BNG + mean:lon_BNG +
                     elev_site + mean:elev_site +
                     (1|siteID), iter=5000, data=CHB)
    b_null <- brm(X_diff~ X_pre +
                      lat_BNG + mean:lat_BNG + lon_BNG + mean:lon_BNG +
                      elev_site + mean:elev_site +
                      (1|siteID), iter=5000, data=CHB)
    
    w_sd <- waic(b_lin_sd)
    w_noSD <- waic(b_lin)
    w_null <- waic(b_null)
    comp <- compare(w_sd, w_noSD)
    comp2 <- compare(w_sd, w_null)
    
    waic_i <- data_frame(species=i, dWAIC_1 = comp[1], dSE_1 = comp[2], dWAIC_2 = comp2[1], dSE_2 = comp2[2])
    predicts_i <- CHB %>% mutate(X_pre = median(X_pre))
    pSD_i <- fitted(b_lin_sd, re_form=NA, newdata=predicts_i) %>% as_data_frame
    pNOSD_i <- fitted(b_lin, re_form=NA, newdata=predicts_i) %>% as_data_frame
    names(pSD_i) <- paste0("sd_", c("est", "SE", "lwr", "upr"))
    names(pNOSD_i) <- paste0("mean_", c("est", "SE", "lwr", "upr"))
    predicts_full_i <- bind_cols(predicts_i, pSD_i, pNOSD_i)
    
    save(b_lin_sd, b_lin, b_null, w_sd, w_noSD, w_null, file = paste0("files/models/", i, "_controls.Rdata"))
    catch_waic[[i]] <- waic_i
    catch_preds[[i]] <- predicts_full_i
    
    saveRDS(catch_waic, "files/waic_updating_controls.rds")
    saveRDS(catch_preds, "files/preds_updating_controls.rds")
}
 