#  --------------------------------------------------------------------------------
# load libraries and setup
# --------------------------------------------------------------------------------

rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(ggplot2)

pop_mc <- data.table::fread(input = here::here("figs/pop_mc.csv"))
pop_mc_means <- data.table::fread(input = here::here("figs/pop_mc_means.csv"))

abm_mc <- data.table::fread(input = here::here("figs/abm_mc.csv"))
abm_mc_means <- data.table::fread(input = here::here("figs/abm_mc_means.csv"))

pop_mc[,"type" := "pop"]
abm_mc[,"type" := "abm"]

mc_dt <- rbind(pop_mc,abm_mc)

abm_mc_means <- melt(abm_mc_means,id.vars = "variable",variable.name = "type")
abm_mc_means[type == "mean", type := "abm"]

pop_mc_means <- melt(pop_mc_means,id.vars = "variable",variable.name = "type")
pop_mc_means[type == "mean", type := "pop"]

mc_means <- merge(abm_mc_means,pop_mc_means,all= TRUE)
setorder(mc_means,"type")

mc_means_dt <- rbind(abm_mc_means,pop_mc_means)

ggplot(data = mc_dt) +
  geom_histogram(aes(value,after_stat(ndensity),fill=variable,linetype=type),position = "identity", color = "black", size = 0.15,alpha=0.6) +
  # geom_vline(data = mc_means,aes(xintercept=value,color=variable,linetype=type)) +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(. ~ variable,scales = "free") +
  theme_bw() +
  theme(axis.title = element_blank())


ggplot(data = merge(mc_dt,mc_means_dt)) +
  geom_violin(aes(x=type,y=value,fill=variable),alpha=0.8) +
  # geom_hline(data = mc_means,aes(yintercept=value,color=variable,linetype=type)) +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(. ~ variable,scales = "free") +
  theme_bw() +
  theme(axis.title = element_blank())
