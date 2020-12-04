# --------------------------------------------------------------------------------
# load libraries and setup
# --------------------------------------------------------------------------------

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
tmax <- 365*100

IC <- calc_equilibrium(NH = NH,X = X)

# dde_odin <- derivs_odin(
#   a = IC$parameters[["a"]],
#   b = IC$parameters[["b"]],
#   c = IC$parameters[["c"]],
#   EIP = IC$parameters[["EIP"]],
#   g = IC$parameters[["g"]],
#   lambda = IC$parameters[["lambda"]],
#   r = IC$parameters[["r"]],
#   LEP = IC$parameters[["LEP"]],
#   IH0 = IC$y0[["IH"]],
#   IV0 = IC$y0[["IV"]],
#   SH0 = IC$y0[["SH"]],
#   SV0 = IC$y0[["SV"]],
#   EH0 = IC$y0[["EH"]],
#   EV0 = IC$y0[["EV"]]
# )
# 
# dde_out <- dde_odin$run(t = 0:tmax)
# 
# dde_out <- as.data.table(dde_out)
# data.table::setnames(x = dde_out,old = "t",new = "time")
# dde_out <- data.table::melt(dde_out,id.vars="time")
# 
# ggplot(data=dde_out) +
#   geom_line(aes(x=time,y=value,color=variable)) +
#   theme_bw()


# --------------------------------------------------------------------------------
# run Gillespie simulation
# --------------------------------------------------------------------------------

Rcpp::sourceCpp(here::here("r-src/exposed/aggregated.cpp"))
Rcpp::sourceCpp(here::here("r-src/discretise.cpp"))

agg_out <- pfsim_aggregated(
  tmax = tmax,
  SH = IC$y0[["SH"]] + IC$y0[["EH"]],
  IH = IC$y0[["IH"]],
  SV = IC$y0[["SV"]] + IC$y0[["EV"]],
  IV = IC$y0[["IV"]],
  parameters = IC$parameters,
  verbose = TRUE
)

agg_out <- as.data.table(discretise(out = agg_out,dt = 1))
agg_out <- data.table::melt(agg_out,id.vars = "time",measure.vars = c("SH","EH","IH","SV","EV","IV"))
agg_out[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
agg_out[, ("sd") := vapply(X = 1:.N,FUN = function(n){sd(value[1:n])},FUN.VALUE = numeric(1)), by = .(variable)]
agg_out[, ("species") := ifelse(variable %in% c("SH","EH","IH"),"Human","Mosquito")]


ic_dt <- data.table::melt(as.data.table(t(as.data.frame(IC$y0))))
colnames(ic_dt)[2] <- "Equilibrium"
ic_dt[, ("species") := ifelse(variable %in% c("SH","EH","IH"),"Human","Mosquito")]

comb_out_agg <- rbind(ic_dt,agg_out,fill=T)

plot_agg <- ggplot(data = comb_out_agg) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.25) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  geom_hline(aes(yintercept=Equilibrium,color=variable),linetype=2,alpha=0.9,size=1.05) +
  facet_wrap(. ~ variable,scales = "free_y") +
  ggtitle("Gillespie simulation vs. DDE") +
  theme_bw()

ggsave(plot = plot_agg, filename = here::here("figs/agg_ts_exact.tiff"),device = "tiff",width = 10,height = 4, compression = "lzw")

plot_agg_sd <- ggplot(data = comb_out_agg) +
  geom_line(aes(x=time,y=mean,color=variable)) + 
  geom_ribbon(aes(x=time,ymin=mean-sd,ymax=mean+sd,fill=variable),alpha=0.25) +
  geom_hline(aes(yintercept=Equilibrium,color=variable),linetype=2,alpha=0.9,size=1.05) +
  facet_wrap(. ~ variable,scales = "free_y") +
  ggtitle("Gillespie simulation vs. DDE") +
  guides(fill=FALSE,color=FALSE) +
  theme_bw()

ggsave(plot = plot_agg_sd, filename = here::here("figs/agg_ts_sd_exact.tiff"),device = "tiff",width = 10,height = 6, compression = "lzw")

plot_agg_hist <- ggplot(data = comb_out_agg) +
  geom_histogram(aes(value,after_stat(density),fill=variable), position = "identity", alpha = 0.5, color = adjustcolor("black",alpha.f = 0.85),size=0.25) + 
  geom_vline(aes(xintercept = Equilibrium, color = variable),linetype=2,alpha=0.9,size=1.05)+
  facet_wrap(. ~ variable,scales = "free") +
  ggtitle("Gillespie simulation vs. DDE") +
  guides(fill=FALSE,color=FALSE) +
  theme_bw()

ggsave(plot = plot_agg_hist, filename = here::here("figs/agg_hist_exact.tiff"),device = "tiff",width = 10,height = 6, compression = "lzw")


# --------------------------------------------------------------------------------
# run MASH simulation
# --------------------------------------------------------------------------------