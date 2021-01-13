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
  geom_histogram(aes(x = value,y = after_stat(density),fill=type),position = "identity", color = "black", size = 0.15,alpha=0.25) +
    facet_wrap(. ~ variable,scales = "free") +
    guides(color = FALSE) +
    theme_bw() +
    theme(axis.title = element_blank())

ggsave(filename = here::here("figs/disagg_hist.tiff"),plot = plot_hist,device = "tiff",height = 6,width = 8,compression = "lzw")

abm[,.(mean = mean(value)),by=variable]
pop[,.(mean = mean(value)),by=variable]