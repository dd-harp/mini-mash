rm(list=ls());gc()
# dev.off()

library(here)
library(data.table)
library(odin)
library(ggplot2)
library(deSolve)
library(Rcpp)

Rcpp::sourceCpp(here::here("src/aggregated/aggregated.cpp"))
Rcpp::sourceCpp(here::here("src/discretise.cpp"))
source(here::here("src/deterministic.R"))

NH <- 1e4
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)

y0_int <- round(IC$y0)

# run DDE model
tmax <- 365*100

dde_odin <- derivs_odin$new(
  a = IC$parameters[["a"]],
  b = IC$parameters[["b"]],
  c = IC$parameters[["c"]],
  EIP = IC$parameters[["EIP"]],
  g = IC$parameters[["g"]],
  lambda = IC$parameters[["lambda"]],
  r = IC$parameters[["r"]],
  LEP = IC$parameters[["LEP"]],
  IH0 = IC$y0[["IH"]],
  IV0 = IC$y0[["IV"]],
  EH0 = IC$y0[["EH"]],
  EV0 = IC$y0[["EV"]],
  SH0 = IC$y0[["SH"]],
  SV0 = IC$y0[["SV"]]
)

dde_out <- dde_odin$run(t = 0:tmax)

dde_out <- as.data.table(dde_out)
data.table::setnames(x = dde_out,old = "t",new = "time")
dde_out <- data.table::melt(dde_out,id.vars = "time",measure.vars = c("SH","EH","IH","SV","EV","IV"))
dde_out[, ("species") := ifelse(variable %in% c("SH","IH"),"Human","Mosquito")]

# plot equilibrium vs DDEs
ggplot(dde_out) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.95) +
  geom_hline(data = eq2dt(IC), mapping = aes(yintercept=value,color=variable),linetype=2,alpha=0.5) +
  facet_wrap(. ~ variable,scales="free_y") +
  theme_bw()

# run aggregated simulation
out <- pfsim_aggregated(
  tmax = tmax,
  SH = y0_int[["SH"]],
  EH = y0_int[["EH"]],
  IH = y0_int[["IH"]],
  SV = y0_int[["SV"]],
  EV = y0_int[["EV"]],
  IV = y0_int[["IV"]],
  parameters = IC$parameters,
  verbose = T
)

out <- data.table::as.data.table(discretise(out = out,dt = 1))
# out <- data.table::as.data.table(out)
out <- data.table::melt(out,id.vars = "time",measure.vars = c("SH","EH","IH","SV","EV","IV"))
out[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out[, ("species") := ifelse(variable %in% c("SH","EH","IH"),"Human","Mosquito")]

plot_agg <- ggplot(data = out) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.5) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.9,data=dde_out,linetype=2) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  facet_wrap(. ~ species,scales = "free_y") +
  ggtitle("Gillespie simulation vs. DDE") +
  theme_bw()

ggsave(plot = plot_agg, filename = here::here("figs/gillespie_dde.tiff"),device = "tiff",width = 14,height = 8)

# run disaggregated population simulation
Rcpp::sourceCpp(here::here("src/disaggregated/disaggregated-pop.cpp"),rebuild = TRUE)

full_out <- run_miniMASH_pop(
  SV = y0_int[["SV"]], EV = y0_int[["EV"]], IV = y0_int[["IV"]], 
  SH = y0_int[["SH"]], EH = y0_int[["EH"]], IH = y0_int[["IH"]],
  parameters = IC$parameters, dt = 5, tmax = tmax 
)
# full_out <- run_miniMASH(parameters = IC$parameters,y0 = round(IC$y0),dt = 5,tmax = tmax)

out_m <- full_out$mosy
colnames(out_m) <- c("time","SV","EV","IV")
out_h <- full_out$human
colnames(out_h) <- c("time","SH","EH","IH")

out_m <- data.table::as.data.table(discretise(out = out_m,dt = 1))
out_m <- data.table::melt(out_m,id.vars="time")
out_m[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out_m[, ("species") := "Mosquito"]

out_h <- data.table::as.data.table(discretise(out = out_h,dt = 1))
out_h <- data.table::melt(out_h,id.vars="time")
out_h[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out_h[, ("species") := "Human"]

plot_disagg <- ggplot(data = rbind(out_h,out_m)) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.5) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.9,data=dde_out,linetype=2) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  facet_wrap(. ~ species,scales = "free_y") +
  ggtitle("MASH vs. DDE") +
  theme_bw()

ggsave(plot = plot_disagg, filename = here::here("figs/MASH_pop_dde.tiff"),device = "tiff",width = 14,height = 8)

# run aggregated population simulation
Rcpp::sourceCpp(here::here("src/disaggregated/disaggregated-abm.cpp"),rebuild = TRUE)

full_out <- run_miniMASH_abm(
  SV = y0_int[["SV"]], EV = y0_int[["EV"]], IV = y0_int[["IV"]], 
  SH = y0_int[["SH"]], EH = y0_int[["EH"]], IH = y0_int[["IH"]],
  parameters = IC$parameters, dt = 5, tmax = tmax 
)

out_m <- full_out$mosy
colnames(out_m) <- c("time","SV","EV","IV")
out_h <- full_out$human
colnames(out_h) <- c("time","SH","EH","IH")

out_m <- data.table::as.data.table(discretise(out = out_m,dt = 1))
out_m <- data.table::melt(out_m,id.vars="time")
out_m[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out_m[, ("species") := "Mosquito"]

out_h <- data.table::as.data.table(discretise(out = out_h,dt = 1))
out_h <- data.table::melt(out_h,id.vars="time")
out_h[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out_h[, ("species") := "Human"]

plot_disagg <- ggplot(data = rbind(out_h,out_m)) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.5) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.9,data=dde_out,linetype=2) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  facet_wrap(. ~ species,scales = "free_y") +
  ggtitle("MASH vs. DDE") +
  theme_bw()

ggsave(plot = plot_disagg, filename = here::here("figs/MASH_abm_dde.tiff"),device = "tiff",width = 14,height = 8)
