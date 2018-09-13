library(dplyr); library(brms)
# tabulate brms output
fnames <- list.files("../mean-variance/files/models/", full.names=T) 
fnames <- fnames[!grepl("controls|JI", fnames)]

bfly_key <- readRDS("files/clean butterfly dataset_v2.rds") %>%
    dplyr::select(species, species_LB) %>% unique

comps_fmtd <- readRDS("files/waic comparisons_formatted.rds") %>% 
    dplyr::select(species, type, dWAIC_1, dSE_1, sp_order, p) %>% 
    unique %>%
    filter(type == "noCont") %>%
    dplyr::select(species, sp_order, dWAIC = dWAIC_1, p)

RMSE <- readRDS("../mean-variance/files/preds_updating.rds") %>%
    bind_rows() %>%
    group_by(species) %>%
    summarise(RMSE = sqrt(mean((X_diff-sd_est)^2)))

sp_list <- gsub(".*models/(.*).Rdata", "\\1", fnames)
sp_list
round_kd <- function(x, digits) formatC(round(x, digits), format='f', digits=digits)
catch_rows <- list()
for(i in seq(fnames)) {
    print(i)
    load(fnames[i])
    sp_i <- sp_list[i]
    summ_i <- summary(b_lin_sd)
    fit_R2 <- bayes_R2(b_lin_sd)
    row_i <- summ_i$fixed %>% 
        .[3:5, c(1,3,4)] %>%
        apply(., c(1,2), round_kd, digits=2) %>%
        as.data.frame %>%
        mutate(est_ul = paste0(Estimate, " (", `l-95% CI`, ", ", `u-95% CI`, ")")) %>%
        dplyr::select(est_ul) %>% 
        t %>%
        as_data_frame %>%
        mutate(species = sp_i, R2 = round_kd(fit_R2[1], 2)) %>%
        select(species, mean = 1, sd = 2, mean_sd = 3, R2) %>%
        left_join(., comps_fmtd) %>%
        left_join(., RMSE) %>%
        dplyr::select(species, sp_order, everything())
    catch_rows[[i]] <- row_i
    
}

model_summary <- bind_rows(catch_rows) %>%
    filter(!species %in% c("Clouded yellow", "Painted lady", "Red admiral")) %>%
    arrange(sp_order) %>%
    mutate(dWAIC = round_kd(dWAIC, 1), 
           RMSE = round_kd(RMSE, 2)) %>% 
    left_join(., bfly_key) %>%
    left_join(., comps_fmtd %>% dplyr::select(species, p)) %>%
    mutate(period = gsub(".*, ", "", sp_order), 
           species_incLB = paste0(1:n(),". ", species, ", ", species_LB))

write.csv(model_summary, "files/model_summary.csv")
