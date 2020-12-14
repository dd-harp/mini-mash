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



tmax <- 1e3
dt <- 5
mc <- 500

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
    Rcpp::sourceCpp(here::here("r-src/exposed/disaggregated-abm.cpp"),showOutput = FALSE)
  }
)

parallel::clusterSetRNGStream(cl = cl,iseed = 32452345L)

system.time(abm_parout <- foreach(i = 1:mc, .export = c("SH","IH","SV","IV","IC","dt","tmax"), .combine = "rbind",.packages = "matrixStats") %dopar% {

  out <- run_miniMASH_abm(SV = SV,IV = IV,SH = SH,IH = IH,parameters = IC$parameters,dt = dt,tmax = tmax)

  # nh <- nrow(out$human)
  # nm <- nrow(out$mosy)
  # ret <- c(
  #   runid = i,
  #   SH = out$human[nh,"SH"],
  #   EH = out$human[nh,"EH"],
  #   IH = out$human[nh,"IH"],
  #   SV = out$mosy[nm,"SV"],
  #   EV = out$mosy[nm,"EV"],
  #   IV = out$mosy[nm,"IV"]
  # )
  
  hmeans <- colWeightedMeans(x = out$human[-nrow(out$human),-1],w = diff(out$human)[,"time"])
  mmeans <- colWeightedMeans(x = out$mosy[-nrow(out$mosy),-1],w = diff(out$mosy)[,"time"])
  ret <- c(i,hmeans,mmeans)
  ret <- setNames(object = ret,nm = c("runid","SH","EH","IH","SV","EV",'IV'))

  return(ret)

})

parallel::stopCluster(cl)


abm_dt <- as.data.table(abm_parout)
abm_dt <- data.table::melt(abm_dt, id.vars = "runid")

ggplot(data=abm_dt) +
  geom_histogram(aes(value,after_stat(density),fill=variable),position = "identity", color = "black", size = 0.1,alpha=0.75) +
  # geom_vline(data = data.frame(variable = names(IC$y0), value = IC$y0),mapping = aes(color = variable,xintercept = value)) +
  facet_wrap(. ~ variable,scales = "free") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.title = element_blank())
