rm(list=ls());gc()

library(here)
library(data.table)
library(ggplot2)
library(deSolve)
library(Rcpp)
library(stocheulerABM)

Rcpp::sourceCpp(here::here("r-src/disaggregated.cpp"))
source(here::here("r-src/deterministic.R"))


NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X,EIP = 11)

tmax <- 365*5

mout <- test_mosquitoes(parameters = IC$parameters,X = X,SV = IC$y0[["SV"]],IV = IC$y0[["IV"]],dt = 5,tmax = tmax)

matplot(mout[,-1],type="l")
IC$parameters[["lambda"]]/IC$parameters[["g"]]
mean(mout[,2])
mean(mout[,3])
mean(mout[,2])+mean(mout[,3])
