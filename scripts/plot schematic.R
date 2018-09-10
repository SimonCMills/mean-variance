# simulate some figures..

## housekeeping ----
library(ggplot2); library(dplyr)
# label info for plot
labs <- data_frame(x = 20, y=8 - .7*(1:4), w_sd = factor(1:4), text = c("sigma == 1", "sigma == 2", "sigma == 3", "sigma == 4"))
labs2 <- data_frame(x = 20, y=8, text = c("sigma == 0"))

## Panel A ----
# fitness function
fitness <- function(x, opt) dnorm(x, opt, 5)*100
# weighted fitness
rf_1 <- function(x, opt=20, w_mean, w_sd) fitness(x, opt)*dnorm(x, w_mean, w_sd)
rf_2 <- function(x, opt=15, w_mean, w_sd) fitness(x, opt)*dnorm(x, w_mean, w_sd)
rf_3 <- function(x, opt=10, w_mean, w_sd) fitness(x, opt)*dnorm(x, w_mean, w_sd)
# realised fitness by integrating across fine-scale responses
intgrt_abbr <- function(x, rf) integrate(rf, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value

# generate responses
wvars <- expand.grid(w_mean = seq(10, 20, length.out=50), w_sd=c(1, 2, 3, 4))
wvars$pred1 <- apply(wvars, 1, intgrt_abbr, rf=rf_1)
wvars$pred2 <- apply(wvars, 1, intgrt_abbr, rf=rf_2)
wvars$pred3 <- apply(wvars, 1, intgrt_abbr, rf=rf_3)

# plot
p1 <- ggplot(wvars, aes(w_mean, pred1, col=factor(w_sd))) + 
    geom_line(aes(y=fitness(w_mean, 20)), col="black", lty=2) +
    geom_line() +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    ylim(0, 8) +
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"))

p2 <- ggplot(wvars, aes(w_mean, pred2, col=factor(w_sd))) + 
    geom_line(aes(y=fitness(w_mean, 15)), col="black", lty=2) +
    geom_line() +
    guides(col=F) +
    scale_colour_brewer(palette="RdBu") +
    theme(legend.position = "bottom") +
    labs(x="Mean temperature", y="Abundance", colour = "SD temperature") +
    ylim(0, 8) +
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"),
          axis.text.y = element_blank(), 
          axis.title.y = element_blank())

p3 <- ggplot(wvars, aes(w_mean, pred3, col=factor(w_sd))) + 
    geom_line(aes(y=fitness(w_mean, 10)), col="black", lty=2) +
    geom_line() +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") +
    ylim(0, 8) +
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    geom_text(data=labs, aes(x, y, label=text), hjust=1, parse=T) +
    geom_text(data=labs2, aes(x, y, label=text), hjust=1, parse=T, inherit.aes=F) +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank())

p_all3 <- egg::ggarrange(p1,p2,p3, nrow=1)
# ggsave("figures/variance_schematic.png",p_all3, width=200, height=80, dpi=300, units="mm")

## Panel B ----
# fitness function
g2_alt <- function(x, mean, sd) {
    p_crit <- 1-pnorm(21, mean, sd)
    p_crit*(dnorm(x, 20, 5)*0) + (1-p_crit)*dnorm(x,20,5)
}
r2_alt <- function(x, w_mean, w_sd) g2_alt(x, mean = w_mean, sd = w_sd)*dnorm(x, mean=w_mean, sd = w_sd)

##
g3_alt <- function(x, mean, sd) {
    p_crit <- 1-pnorm(21, mean, sd)
    p_crit*(dnorm(x, 15, 5)*0) + (1-p_crit)*dnorm(x,15,5)
}
r3_alt <- function(x, w_mean, w_sd) g3_alt(x, mean = w_mean, sd = w_sd)*dnorm(x, mean=w_mean, sd = w_sd)

g4_alt <- function(x, mean, sd) {
    p_crit <- 1-pnorm(21, mean, sd)
    p_crit*(dnorm(x, 10, 5)*0) + (1-p_crit)*dnorm(x,10,5)
}
r4_alt <- function(x, w_mean, w_sd) g4_alt(x, mean = w_mean, sd = w_sd)*dnorm(x, mean=w_mean, sd = w_sd)

wvars$p2_alt <- apply(wvars, 1, function(x) integrate(r2_alt, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p3_alt <- apply(wvars, 1, function(x) integrate(r3_alt, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p4_alt <- apply(wvars, 1, function(x) integrate(r4_alt, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)

p1_alt <- ggplot(wvars, aes(w_mean, p2_alt*100, col=factor(w_sd))) + 
    geom_line(aes(y=fitness(w_mean, 20)), col="black", lty=2) +
    geom_line()  +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    ylim(0, 8)+
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme(panel.background = element_rect(colour="black")) 

##
p2_alt <- ggplot(wvars, aes(w_mean, p3_alt*100, col=factor(w_sd))) + 
    geom_line(aes(y=fitness(w_mean, 15)), col="black", lty=2) +
    geom_line()  +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    ylim(0, 8)+
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank())

##
p3_alt <- ggplot(wvars, aes(w_mean, p4_alt*100, col=factor(w_sd))) + 
    geom_line(aes(y=fitness(w_mean, 10)), col="black", lty=2) +
    geom_line()  +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    ylim(0, 8)+
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank())
p_all3 <- egg::ggarrange(p1_alt,p2_alt,p3_alt, nrow=1)
# ggsave("figures/variance_schematic_alt.png",p_all3, width=200, height=80, dpi=300, units="mm")

p1_labelled <- p1 + labs(title="(a) fixed response") +
    theme(plot.title = element_text(hjust=0, size=10, face="bold"))

p1_alt_labelled <- p1_alt + labs(title="(b) impaired response when temperature threshold exceeded") +
    theme(plot.title = element_text(hjust=0, size=10, face="bold"))

p_all3 <- egg::ggarrange(p1_labelled, p2, p3, p1_alt_labelled, p2_alt,p3_alt, nrow=2)
ggsave("figures/variance_schematic_all.png",p_all3, width=200, height=130, dpi=300, units="mm")