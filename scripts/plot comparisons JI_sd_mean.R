# plot 'empirical' models
library(dplyr); library(ggplot2)

# read in model comparison waic scores
fmtd <- readRDS("../mean-variance/files/waic comparisons_formatted.rds") %>%
    filter(type=="cont")
comps <- readRDS("../mean-variance/files/comps_JI_brms_updating_v2.rds") %>%
    bind_rows()
comps_full <- comps %>%
    left_join(., fmtd %>% dplyr::select(species, sp_order)) %>%
    mutate(p_JI_mean = (1-pnorm(0, dWAIC_JI_mean, SE_JI_mean)) < .05, 
           p_sd_mean = (1-pnorm(0, dWAIC_sd_mean, SE_sd_mean)) < .05, 
           p_class = case_when(p_JI_mean & !p_sd_mean~ "empirical only", 
                               !p_JI_mean & p_sd_mean~ "sd only", 
                               p_JI_mean & p_sd_mean~ "both", 
                               TRUE ~ "neither"), 
           p_class = factor(p_class, levels=c("neither", "sd only", "both", "empirical only"))) 

# rep plot for latter models ----
p2 <- ggplot(comps_full, aes(y=dWAIC_JI_mean, 
                             x=sp_order, 
                             ymin=dWAIC_JI_mean- 1.96*SE_JI_mean, 
                             ymax = dWAIC_JI_mean + 1.96*SE_JI_mean)) +
    geom_errorbar(width=0) +
    geom_point(pch=21, size=2, fill="black") +
    coord_flip() +
    geom_vline(xintercept=0, lty=2) +

    labs(y = bquote(Delta*WAIC), x="") +
    geom_hline(yintercept=0, lty=2) +
    theme_classic() +
    theme(axis.text= element_text(colour="black"),
          panel.background = element_rect(colour="black"),
          axis.title.y = element_blank()) +
    guides(col=F, fill=F)

ggsave("figures/modComp_JI_vsNull.png",plot=p2, width=120, height=150, units="mm", dpi=200)

# 2d plot ----
p3 <- ggplot(comps_full, aes(dWAIC_sd_mean, dWAIC_JI_mean, fill=p_class)) + 
    scale_fill_manual(values=c("black", RColorBrewer::brewer.pal(9, name = "RdBu")[c(2,8)])) + 
    geom_vline(xintercept=0, lty=2) +
    geom_hline(yintercept=0, lty=2) +
    geom_abline() +
    geom_point(pch=21, size=2) + 
    xlim(c(-75, 5)) + ylim(c(-75, 5)) +
    coord_equal() +
    labs(x = bquote(Delta*WAIC[Phen]), y=bquote(Delta*WAIC["Emp"])) +
    theme_classic() +
    guides(fill=F) +
    theme(panel.background = element_rect(colour="black"))

ggsave("figures/modComp_JI_vsPhen_2d.png", p3, width=100, height=100, units="mm", dpi=200)