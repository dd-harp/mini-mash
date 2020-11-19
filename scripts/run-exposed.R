rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(odin)
library(ggplot2)
library(deSolve)
library(Rcpp)

Rcpp::sourceCpp(here::here("r-src/exposed/aggregated.cpp"))
Rcpp::sourceCpp(here::here("r-src/discretise.cpp"))
source(here::here("r-src/deterministic.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)

tmax <- 365*100

dde_odin <- derivs_odin(
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
  SH0 = IC$y0[["SH"]],
  SV0 = IC$y0[["SV"]]
)

dde_out <- dde_odin$run(t = 0:tmax)
# dde_out <- dede(y = IC$y0,times = 0:tmax,func = derivs,parms = IC$parameters)

dde_out <- as.data.table(dde_out)
data.table::setnames(x = dde_out,old = "t",new = "time")
dde_out <- data.table::melt(dde_out,id.vars = "time",measure.vars = c("SH","IH","SV","IV"))
dde_out[, ("species") := ifelse(variable %in% c("SH","IH"),"Human","Mosquito")]

out <- pfsim_aggregated(
  tmax = tmax,
  SH = IC$y0[["SH"]],
  IH = IC$y0[["IH"]],
  SV = IC$y0[["SV"]],
  IV = IC$y0[["IV"]],
  parameters = IC$parameters,
  verbose = T
)

out <- data.table::as.data.table(discretise(out = out,dt = 1))
out <- data.table::melt(out,id.vars = "time",measure.vars = c("SH","EH","IH","SV","EV","IV"))
out[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out[, ("species") := ifelse(variable %in% c("SH","EH","IH"),"Human","Mosquito")]

ggplot(data = out) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.5) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.9,linetype="dashed",size=1.15,data=dde_out) +
  facet_wrap(. ~ species,scales = "free_y") +
  ggtitle("Gillespie simulation vs. DDE") +
  theme_bw()
