rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(ggplot2)
library(deSolve)
library(Rcpp)
library(stocheulerABM)

Rcpp::sourceCpp(here::here("r-src/disaggregated.cpp"))
Rcpp::sourceCpp(here::here("r-src/disaggregateddebug.cpp"))
Rcpp::sourceCpp(here::here("r-src/discretise.cpp"))
source(here::here("r-src/deterministic.R"))


NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X,EIP = 11)
IC$y0

tmax <- 365*100

mout <- test_mosquitoes(parameters = IC$parameters,X = X,SV = IC$y0[["SV"]],IV = IC$y0[["IV"]],dt = 5,tmax = tmax)
mout <- as.data.table(discretise(out = mout,dt = 1))
mout <- data.table::melt(data = mout,id.vars="time")

ggplot(data = mout) +
  geom_line(aes(x=time,y=value,color=variable)) +
  geom_hline(yintercept = IC$y0[3:4],linetype=2,alpha=0.75) +
  theme_bw()

mout[variable == "S", mean(value)]
mout[variable == "I", mean(value)]


hout <- test_humans(parameters = IC$parameters,SH = IC$y0[["SH"]],IH = IC$y0[["IH"]],IV = IC$y0[["IV"]],dt = 5,tmax = tmax)
hout <- as.data.table(discretise(out = hout,dt = 1))
hout <- data.table::melt(data = hout,id.vars="time")

ggplot(data = hout) +
  geom_line(aes(x=time,y=value,color=variable)) +
  geom_hline(yintercept = IC$y0[1:2],linetype=2,alpha=0.75) +
  theme_bw()

hout[variable == "S", mean(value)]
hout[variable == "I", mean(value)]
