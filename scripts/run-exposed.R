rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(odin)
library(ggplot2)
library(deSolve)
library(Rcpp)

Rcpp::sourceCpp(here::here("r-src/exposed/aggregated.cpp"))
Rcpp::sourceCpp(here::here("r-src/discretise.cpp"))
source(here::here("r-src/deterministic.R"))

NH <- 1e3
X <- 0.3

IC <- calc_equilibrium(NH = NH,X = X)

tmax <- 365*100

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
#   SV0 = IC$y0[["SV"]]
# )
# 
# dde_out <- dde_odin$run(t = 0:tmax)
# 
# dde_out <- as.data.table(dde_out)
# data.table::setnames(x = dde_out,old = "t",new = "time")
# dde_out <- data.table::melt(dde_out,id.vars = "time",measure.vars = c("SH","IH","SV","IV"))
# dde_out[, ("species") := ifelse(variable %in% c("SH","IH"),"Human","Mosquito")]

ic_dt <- data.table::melt(as.data.table(t(as.data.frame(IC$y0))))
colnames(ic_dt)[2] <- "Equilibrium"
ic_dt[, ("species") := ifelse(variable %in% c("SH","EH","IH"),"Human","Mosquito")]

out <- pfsim_aggregated(
  tmax = tmax,
  SH = IC$y0[["SH"]],
  IH = IC$y0[["IH"]],
  SV = IC$y0[["SV"]],
  IV = IC$y0[["IV"]],
  parameters = IC$parameters,
  verbose = T
)

out <- data.table::as.data.table(discretise(out = out,dt = 1))
out <- data.table::melt(out,id.vars = "time",measure.vars = c("SH","EH","IH","SV","EV","IV"))
out[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out[, ("species") := ifelse(variable %in% c("SH","EH","IH"),"Human","Mosquito")]

# plot_agg <- ggplot(data = rbind(ic_dt,out,fill=T)) +
#   geom_line(aes(x=time,y=value,color=variable),alpha=0.25) +
#   geom_line(aes(x=time,y=mean,color=variable)) +
#   geom_hline(aes(yintercept=Equilibrium,color=variable),linetype=2,alpha=0.9,size=1.05) +
#   facet_wrap(. ~ species,scales = "free_y") +
#   ggtitle("Gillespie simulation vs. DDE") +
#   theme_bw()

comb_out_agg <- rbind(ic_dt,out,fill=T)

plot_agg <- ggplot(data = comb_out_agg[!variable %in% c("EV","EH"), ]) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.25) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  geom_hline(aes(yintercept=Equilibrium,color=variable),linetype=2,alpha=0.9,size=1.05) +
  facet_wrap(. ~ variable,scales = "free_y") +
  ggtitle("Gillespie simulation vs. DDE") +
  theme_bw()

ggsave(plot = plot_agg, filename = here::here("figs/gillespie_dde_E.tiff"),device = "tiff",width = 12,height = 6)

Rcpp::sourceCpp(here::here("r-src/exposed/disaggregated.cpp"))

full_out <- run_miniMASH(parameters = IC$parameters,y0 = round(IC$y0),dt = 5,tmax = tmax)

out_m <- data.table::as.data.table(discretise(out = full_out$mosquito,dt = 1))
out_m <- data.table::melt(out_m,id.vars="time")
out_m[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out_m[, ("species") := "Mosquito"]

out_h <- data.table::as.data.table(discretise(out = full_out$human,dt = 1))
out_h <- data.table::melt(out_h,id.vars="time")
out_h[, ("mean") := cumsum(value)/1:.N, by = .(variable)]
out_h[, ("species") := "Human"]

# plot_disagg <- ggplot(data = rbind(ic_dt,rbind(out_h,out_m),fill=T)) +
#   geom_line(aes(x=time,y=value,color=variable),alpha=0.25) +
#   geom_line(aes(x=time,y=mean,color=variable)) +
#   geom_hline(aes(yintercept=Equilibrium,color=variable),linetype=2,alpha=0.9,size=1.05) +
#   facet_wrap(. ~ species,scales = "free_y") +
#   ggtitle("MASH vs. DDE") +
#   theme_bw()

comb_out_disagg <- rbind(ic_dt,rbind(out_h,out_m),fill=T)

plot_disagg <- ggplot(data = comb_out_disagg[!variable %in% c("EV","EH"),]) +
  geom_line(aes(x=time,y=value,color=variable),alpha=0.25) +
  geom_line(aes(x=time,y=mean,color=variable)) +
  geom_hline(aes(yintercept=Equilibrium,color=variable),linetype=2,alpha=0.9,size=1.05) +
  facet_wrap(. ~ variable,scales = "free_y") +
  ggtitle("MASH vs. DDE") +
  theme_bw()

ggsave(plot = plot_disagg, filename = here::here("figs/MASH_dde_E.tiff"),device = "tiff",width = 12,height = 6)


# plot histograms

out <- pfsim_aggregated(
  tmax = tmax,
  SH = IC$y0[["SH"]],
  IH = IC$y0[["IH"]],
  SV = IC$y0[["SV"]],
  IV = IC$y0[["IV"]],
  parameters = IC$parameters,
  verbose = T
)

out <- melt(as.data.table(out),id.vars="time")
ggplot(data = out) +
  geom_histogram(aes(value,after_stat(density),fill=variable,color=variable,alpha=0.75)) +
  facet_wrap(. ~ variable,scales="free")+
  guides(alpha=FALSE)+
  theme_bw()
