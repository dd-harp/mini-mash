#  --------------------------------------------------------------------------------
# load libraries and setup
# --------------------------------------------------------------------------------

rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(ggplot2)
library(deSolve)
library(Rcpp)
library(matrixStats)

library(foreach)
library(doParallel)

source(here::here("r-src/deterministic-exposed.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)

tmax <- 5e4
dt <- 5

abm_IC <- round(IC$y0)
SH <- abm_IC[["SH"]]
IH <- sum(abm_IC[c("EH","IH")])
SV <- abm_IC[["SV"]]
IV <- sum(abm_IC[c("EV","IV")])

Rcpp::sourceCpp(here::here("r-src/exposed/disaggregated-pop.cpp"))

system.time(
  simout <- run_miniMASH_pop(SV = SV,IV = IV,SH = SH,IH = IH,parameters = IC$parameters,dt = dt,tmax = tmax)
)

hout <- simout$human
hout <- hout[hout[,"time"] > 500, ]

hmeans <- colWeightedMeans(x = hout[-nrow(hout),2:4],w = diff(hout[,"time"]))
IC$y0[1:3]
