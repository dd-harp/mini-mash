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

source(here::here("r-src/deterministic-exposed.R"))

NH <- 1e3
X <- 0.3
tmax <- 365*100

IC <- calc_equilibrium(NH = NH,X = X)


Rcpp::sourceCpp(here::here("r-src/exposed/test-human-abm.cpp"),showOutput = FALSE)


hum_IC <- round(IC$y0[1:3] )

hout <- test_humans(SH = hum_IC[1],IH = sum(hum_IC[2:3]),IV = IC$y0[["IV"]],parameters = IC$parameters,dt = 5,tmax = 1e3)

state_df <- as.data.frame(hout)
inf_per <- NULL
exp_per <- NULL
for(i in 1:(nrow(state_df)-1)){
  stateseq <- state_df[i:(i+1),"states"]
  if(all(stateseq == c("E","I"))){
    exp_per <- c(exp_per, (state_df[(i+1),"times"] - state_df[i,"times"]) )
  }
  if(all(stateseq == c("I","S"))){
    inf_per <- c(inf_per, (state_df[(i+1),"times"] - state_df[i,"times"]) )
  }
}

hist(inf_per[inf_per>0])
mean(inf_per[inf_per>0])

state_dt <- as.data.table(state_df)
