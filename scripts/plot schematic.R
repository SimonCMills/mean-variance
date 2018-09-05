# simulate some figures..
library(ggplot2)
# g1 <- function(x) exp(1 + .1*x)
# r1 <- function(x, w_mean, w_sd) g1(x)*dnorm(x, mean=w_mean, sd = w_sd)

g2 <- function(x) dnorm(x, 20, 5)
r2 <- function(x, w_mean, w_sd) g2(x)*dnorm(x, mean=w_mean, sd = w_sd)

g3 <- function(x) dnorm(x, 15, 5)
g2_alt <- function(x, w_mean, w_sd) ifelse(pnorm(10, w_mean, w_sd) < .01, dnorm(x, 20, 5), dnorm(x, 20, 5))
r2_alt <- function(x, w_mean, w_sd) g2_alt(x, w_mean, w_sd)*dnorm(x, mean=w_mean, sd = w_sd)
r3 <- function(x, w_mean, w_sd) g3(x)*dnorm(x, mean=w_mean, sd = w_sd)

g4 <- function(x) dnorm(x, 10, 5)
r4 <- function(x, w_mean, w_sd) g4(x)*dnorm(x, mean=w_mean, sd = w_sd)

x <- seq(-10, 10, .5)
plot(r2(x, 10, 5))
qnorm(.2, 0, 1)


wvars <- expand.grid(w_mean = seq(10, 20, length.out=50), w_sd=c(1, 2, 3, 4))
# wvars$p2 <- apply(wvars, 1, function(x) integrate(r2, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p2 <- apply(wvars, 1, function(x) integrate(r2, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p3 <- apply(wvars, 1, function(x) integrate(r3, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p4 <- apply(wvars, 1, function(x) integrate(r4, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)

g2_alt <- function(x, w_mean, w_sd) ifelse(pnorm(10, w_mean, w_sd) < .01, dnorm(x, 20, 5), dnorm(x, 20, 5))
r2_alt <- function(x, w_mean, w_sd) g2_alt(x, w_mean, w_sd)*dnorm(x, mean=w_mean, sd = w_sd)
wvars$p2_alt <- apply(wvars, 1, function(x) integrate(r2_alt, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p3 <- apply(wvars, 1, function(x) integrate(r3, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p4 <- apply(wvars, 1, function(x) integrate(r4, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)




# wvars$p2 <- apply(wvars, 1, function(x) quantile(g2(rnorm(1e5, mean=x[1], sd=x[2])), .8))
# wvars$p2_2 <- apply(wvars, 1, function(x) quantile(g2(rnorm(1e5, mean=x[1], sd=x[2])), .2))
# wvars$p3 <- apply(wvars, 1, function(x) integrate(r3, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
# wvars$p4 <- apply(wvars, 1, function(x) quantile(g4(rnorm(1e5, mean=x[1], sd=x[2])), .8))

p1 <- ggplot(wvars, aes(w_mean, p2*100, col=factor(w_sd))) + 
    geom_line(aes(y=g2(w_mean)*100), col="black", lty=2) +
    geom_line() +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    ylim(0, 8)+
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"))

p2 <- ggplot(wvars, aes(w_mean, p3*100, col=factor(w_sd))) + 
    geom_line(aes(y=100*g3(w_mean)), col="black", lty=2) +
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

labs <- data_frame(x = 20, y=8 - .7*(1:4), w_sd = factor(1:4), text = c("sigma == 1", "sigma == 2", "sigma == 3", "sigma == 4"))
labs2 <- data_frame(x = 20, y=8, text = c("sigma == 0"))
p3 <- ggplot(wvars, aes(w_mean, p4*100, col=factor(w_sd))) + 
    geom_line(aes(y=100*g4(w_mean)), col="black", lty=2) +
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
ggsave("figures/variance_schematic.png",p_all3, width=200, height=80, dpi=300, units="mm")


##
g2_alt <- function(x, mean, sd) {
    p_crit <- 1-pnorm(21, mean, sd)
    p_crit*(dnorm(x, 20, 5)*.1) + (1-p_crit)*dnorm(x,20,5)
}
r2_alt <- function(x, w_mean, w_sd) g2_alt(x, mean = w_mean, sd = w_sd)*dnorm(x, mean=w_mean, sd = w_sd)

##
g3_alt <- function(x, mean, sd) {
    p_crit <- 1-pnorm(21, mean, sd)
    p_crit*(dnorm(x, 15, 5)*.1) + (1-p_crit)*dnorm(x,15,5)
}
r3_alt <- function(x, w_mean, w_sd) g3_alt(x, mean = w_mean, sd = w_sd)*dnorm(x, mean=w_mean, sd = w_sd)

g4_alt <- function(x, mean, sd) {
    p_crit <- 1-pnorm(21, mean, sd)
    TPC <- 1000*dnorm(x, 10, 5)
    p_crit*(TPC*.1) + (1-p_crit)*TPC
}
r4_alt <- function(x, w_mean, w_sd) g4_alt(x, mean = w_mean, sd = w_sd)*dnorm(x, mean=w_mean, sd = w_sd)

wvars$p2_alt <- apply(wvars, 1, function(x) integrate(r2_alt, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p3_alt <- apply(wvars, 1, function(x) integrate(r3_alt, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)
wvars$p4_alt <- apply(wvars, 1, function(x) integrate(r4_alt, -100, 100, w_mean=x[1], w_sd=x[2])[1]$value)

1-pnorm(21, 21, 1)

p1_alt <- ggplot(wvars, aes(w_mean, p2_alt*100, col=factor(w_sd))) + 
    geom_line(aes(y=g2(w_mean)*100), col="black", lty=2) +
    geom_line()  +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    # ylim(0, 8)+
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme(panel.background = element_rect(colour="black")) 

##
p2_alt <- ggplot(wvars, aes(w_mean, p3_alt*100, col=factor(w_sd))) + 
    geom_line(aes(y=g3(w_mean)*100), col="black", lty=2) +
    geom_line()  +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    # ylim(0, 8)+
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank())

##
p3_alt <- ggplot(wvars, aes(w_mean, p4_alt, col=factor(w_sd))) + 
    geom_line(aes(y=g4(w_mean)*1000), col="black", lty=2) +
    geom_line(aes(y=g4(w_mean)*1000*.1), col="black", lty=2) +
    geom_line()  +
    scale_colour_brewer(palette="RdBu") +
    guides(col=F) +
    labs(x="Mean temperature", y="Abundance") + 
    # ylim(0, 8)+
    scale_x_continuous(breaks=seq(0, 30, 2)) +
    theme_classic() +
    theme_classic() +
    theme(panel.background = element_rect(colour="black"), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank())
p3_alt
p_all3 <- egg::ggarrange(p1_alt,p2_alt,p3_alt, nrow=1)
ggsave("figures/variance_schematic_alt.png",p_all3, width=200, height=80, dpi=300, units="mm")

p1_labelled <- p1 + labs(title="(a) fixed response") +
    theme(plot.title = element_text(hjust=0, size=10, face="bold"))

p1_alt_labelled <- p1_alt + labs(title="(b) impaired response when temperature threshold exceeded") +
    theme(plot.title = element_text(hjust=0, size=10, face="bold"))


p_all3 <- egg::ggarrange(p1_labelled, p2, p3, p1_alt_labelled, p2_alt,p3_alt, nrow=2)
ggsave("figures/variance_schematic_all.png",p_all3, width=200, height=130, dpi=300, units="mm")