rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(ggplot2)
library(deSolve)
library(Rcpp)
library(stocheulerABM)

source(here::here("r-src/deterministic-exposed.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)



tmax <- 50

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
  EH0 = 250,
  EV0 = 0,
  SH0 = IC$y0[["SH"]],
  SV0 = IC$y0[["SV"]]
)

dde_out <- dde_odin$run(t = 0:tmax)

dde_out <- as.data.table(dde_out)
data.table::setnames(x = dde_out,old = "t",new = "time")
dde_out <- data.table::melt(dde_out,id.vars = "time")

ggplot(data = dde_out) +
  geom_line(aes(x=time,y=value,color=variable),size=1.05) +
  theme_bw()

dde_out[variable=="EH",]
dde_out[variable=="EV",]
