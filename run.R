rm(list=ls());gc()

library(here)
library(data.table)
library(ggplot2)
library(deSolve)
library(Rcpp)
library(stocheulerABM)

Rcpp::sourceCpp(here::here("r-src/aggregated.cpp"))
source(here::here("r-src/deterministic.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)


tmax <- 365*50

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
out <- data.table::melt(out,id.vars = "time",measure.vars = c("SH","IH","SV","IV"))
out[, ("mean") := cumsum(value)/1:.N, by = .(variable)]

# ggplot(data = out) +
#   geom_line(aes(x=time,y=value,color=variable)) +
#   theme_bw()

dde_out <- dede(y = IC$y0,times = 0:tmax,func = derivs,parms = IC$parameters)

dde_out <- as.data.table(dde_out)
dde_out <- data.table::melt(dde_out,id.vars = "time",measure.vars = c("SH","IH","SV","IV"))

ggplot(data = out) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.75) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.75,data=dde_out,linetype=2) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  theme_bw()

