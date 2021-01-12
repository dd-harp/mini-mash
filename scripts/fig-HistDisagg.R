# --------------------------------------------------------------------------------
#
#   Compare disaggregated-pop and disaggregated-abm output
#   Sean L. Wu (slwood89@gmail.com)
#   January 2021
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(ggplot2)
library(odin)

# DDE solutions
source(here::here("src/deterministic.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)

abm <- data.table::fread(file = here::here("data/abm_means.csv"))
pop <- data.table::fread(file = here::here("data/pop_means.csv"))

abm[,"type" := "abm"]
pop[,"type" := "pop"]

plot_hist <- ggplot(data = rbind(abm,pop)) +
  geom_histogram(aes(x = value,y = after_stat(density),fill=interaction(variable,type)),position = "identity", color = "black", size = 0.15,alpha=0.) +
    facet_wrap(. ~ variable,scales = "free") +
    guides(fill = FALSE, color = FALSE) +
    theme_bw() +
    theme(axis.title = element_blank())