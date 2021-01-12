# --------------------------------------------------------------------------------
#
#   Run disaggregated-abm and compare to DDE solutions
#   Sean L. Wu (slwood89@gmail.com)
#   January 2021
#
# --------------------------------------------------------------------------------

#  --------------------------------------------------------------------------------
# load libraries and setup
# --------------------------------------------------------------------------------

rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(ggplot2)
library(odin)
library(Rcpp)
library(matrixStats)

library(foreach)
library(doParallel)

# DDE solutions
source(here::here("src/deterministic.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)

# simulation parameters
tmax <- 2e3
dt <- 5
mc <- 5e3

# initial conditions
abm_IC <- round(IC$y0)
SH <- abm_IC[["SH"]]
IH <- sum(abm_IC[c("EH","IH")])
SV <- abm_IC[["SV"]]
IV <- sum(abm_IC[c("EV","IV")])

cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

parallel::clusterEvalQ(
  cl = cl,
  expr = {
    Rcpp::sourceCpp(here::here("src/disaggregated/disaggregated-abm.cpp"),showOutput = FALSE)
  }
)

parallel::clusterSetRNGStream(cl = cl,iseed = 32452345L)

system.time(abm_parout <- foreach(i = 1:mc, .export = c("SH","IH","SV","IV","IC","dt","tmax"), .combine = "rbind",.packages = "matrixStats") %dopar% {

  out <- run_miniMASH_abm(SV = SV,IV = IV,SH = SH,IH = IH,parameters = IC$parameters,dt = dt,tmax = tmax)

  out$human <- out$human[out$human[,"time"] > 500,]
  out$mosy <- out$mosy[out$mosy[,"time"] > 500,]
  hmeans <- colWeightedMeans(x = out$human[-nrow(out$human),-1],w = diff(out$human)[,"time"])
  mmeans <- colWeightedMeans(x = out$mosy[-nrow(out$mosy),-1],w = diff(out$mosy)[,"time"])
  ret <- c(i,hmeans,mmeans)
  ret <- setNames(object = ret,nm = c("runid","SH","EH","IH","SV","EV",'IV'))

  return(ret)

})

parallel::stopCluster(cl)

abm_dt <- as.data.table(abm_parout)
abm_dt <- data.table::melt(abm_dt, id.vars = "runid")

data.table::fwrite(x = abm_dt,file = here::here("data/abm_means.csv"))
