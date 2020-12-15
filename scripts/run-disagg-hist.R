#  --------------------------------------------------------------------------------
# load libraries and setup
# ---------------------------------------------------------------------------------

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

tmax <- 1e4
dt <- 5
mc <- 5e3

abm_IC <- round(IC$y0)
SH <- abm_IC[["SH"]]
IH <- sum(abm_IC[c("EH","IH")])
SV <- abm_IC[["SV"]]
IV <- sum(abm_IC[c("EV","IV")])

y0 <- c("SH"=SH,"IH"=IH,"SV"=SV,"IV"=IV)


#  --------------------------------------------------------------------------------
# inexact
# --------------------------------------------------------------------------------

cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

parallel::clusterEvalQ(
  cl = cl,
  expr = {
    Rcpp::sourceCpp(here::here("r-src/exposed/disaggregated.cpp"),showOutput = FALSE)
  }
)

parallel::clusterSetRNGStream(cl = cl,iseed = 96841588L)

system.time(inexact_parout <- foreach(i = 1:mc, .export = c("SH","IH","SV","IV","IC","dt","tmax","y0"), .combine = "rbind",.packages = "matrixStats") %dopar% {

  out <- run_miniMASH_exactbm(y0 = y0,parameters = IC$parameters,dt = dt,tmax = tmax,verbose = FALSE)

  out$human <- out$human[out$human[,"time"] > 500,]
  out$mosy <- out$mosy[out$mosy[,"time"] > 500,]
  hmeans <- colMeans(out$human[,2:4])
  mmeans <- colMeans(out$mosquito[,2:4])
  ret <- c(i,hmeans,mmeans)
  ret <- setNames(object = ret,nm = c("runid","SH","EH","IH","SV","EV",'IV'))

  return(ret)

})

parallel::stopCluster(cl)


inexact_dt <- as.data.table(inexact_parout)
inexact_dt <- data.table::melt(inexact_dt, id.vars = "runid")

data.table::fwrite(x = inexact_dt,file = here::here("figs/pop_mc.csv"))

inexact_dt_sum <- inexact_dt[,.(mean = mean(value)),by = variable]
inexact_dt_sum$dde <- IC$y0

data.table::fwrite(x = inexact_dt_sum,file = here::here("figs/pop_mc_means.csv"))

plot_inexact_hist <- ggplot(data=merge(inexact_dt,inexact_dt_sum)) +
  geom_histogram(aes(value,after_stat(density),fill=variable),position = "identity", color = "black", size = 0.15,alpha=0.6) +
  geom_vline(aes(color = variable,xintercept = mean),linetype=3,size=0.5) +
  geom_vline(aes(color = variable,xintercept = dde),linetype=2,size=0.8) +
  facet_wrap(. ~ variable,scales = "free") +
  guides(fill = FALSE, color = FALSE) +
  theme_bw() +
  theme(axis.title = element_blank())

ggsave(plot = plot_inexact_hist, filename = here::here("figs/disaggregated_hist.tiff"),device = "tiff",width = 10,height = 6, compression = "lzw")