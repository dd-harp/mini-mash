rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(odin)
library(ggplot2)
library(Rcpp)
library(parallel)
library(foreach)
library(doParallel)

source(here::here("r-src/deterministic.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)

tmax <- 365*25
nrep <- 1e4

# set up cluster and source the file on each core
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 50694091L)

parallel::clusterEvalQ(cl,{
  Rcpp::sourceCpp(here::here("r-src/aggregated.cpp"))
  Rcpp::sourceCpp(here::here("r-src/disaggregated-set.cpp"))
})

system.time(finalstate <- foreach(i = 1:nrep,.combine = "rbind",.export = c("IC","tmax")) %dopar% {

  outG <- pfsim_aggregated(
    tmax = tmax,
    SH = IC$y0[["SH"]],
    IH = IC$y0[["IH"]],
    SV = IC$y0[["SV"]],
    IV = IC$y0[["IV"]],
    parameters = IC$parameters,
    verbose = FALSE
  )

  outM <- run_miniMASH(
    parameters = IC$parameters,
    y0 = round(IC$y0),
    dt = 5,
    tmax = tmax
  )

  outG <- tail(outG,1)[,c("SH","IH","SV","IV")]
  outM <- c(tail(outM$human,1)[,c("S","I")],tail(outM$mosquito,1)[,c("S","I")])
  outM <- setNames(outM,c("SH","IH","SV","IV"))
  out <- as.data.frame(rbind(outG,outM))
  out$type <- c("Gillespie","MASH")
  return(out)

})

# clean up the parallel cluster and remove it
parallel::stopCluster(cl);rm(cl);gc()

finalstatedt <- as.data.table(finalstate)
finalstatedt <- melt(finalstatedt,id.vars = "type")

plot_hist <- ggplot(data =finalstatedt) +
  geom_histogram(aes(x=value,y=after_stat(density),fill=type),alpha=0.5,bins=50,position="identity") +
  facet_wrap(. ~ variable,scales = "free") +
  guides(fill = FALSE) +
  ggtitle("State Variable Comparison") +
  theme_bw() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), strip.text.x = element_text(size = 14))

ggsave(plot = plot_hist, filename = here::here("figs/hist_compare.tiff"),device = "tiff",width = 10,height = 8)
