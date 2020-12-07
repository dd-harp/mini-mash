# --------------------------------------------------------------------------------
# load libraries and setup
# --------------------------------------------------------------------------------

rm(list=ls());gc()
dev.off()
library(ggplot2)


TWICE <- 5
nborn <- 5
lambda <- nborn/TWICE
g <- 1/10


S0 <- 6

S0_mosy <- data.frame(
  enter = rep(0,S0),
  leave = pmin(rexp(n = S0,rate = g),TWICE)
)

emerge <- rpois(n = 1,lambda = lambda*TWICE)
emerge_t <- runif(n = emerge,min = 0,max = TWICE)

emerge_mosy <- data.frame(
  enter = emerge_t,
  leave = pmin(emerge_t + rexp(n = emerge,rate = g),TWICE)
)

# individual histories
indiv_traces <- rbind(S0_mosy,emerge_mosy)
indiv_traces <- indiv_traces[order(indiv_traces$leave),]
indiv_traces$id <- as.integer(1:nrow(indiv_traces))

# aggregated ensemble
ensemble_trace <- data.frame(
  time = 0,count = sum(indiv_traces$enter < 2e-16),event = "begin"
)

trace_temp <- data.frame(
  time = c(indiv_traces$enter,indiv_traces$leave),
  delta = c(rep(1,times=length(indiv_traces$enter)),rep(-1,times=length(indiv_traces$enter))),
  event = c(paste0(indiv_traces$id," emerges"),paste0(indiv_traces$id," dies"))
)
trace_temp <- trace_temp[-which(trace_temp$time==0.0),]
trace_temp <- trace_temp[-which(trace_temp$time==5.0),]
trace_temp <- trace_temp[order(trace_temp$time),]

while(nrow(trace_temp) > 0){
  ensemble_trace <- rbind(
    ensemble_trace,
    data.frame(time = trace_temp[1,'time'], count = ensemble_trace[nrow(ensemble_trace),"count"] + trace_temp[1,"delta"], event = trace_temp[1,"event"])
  )
  trace_temp <- trace_temp[-1,]
}

ensemble_trace <- rbind(
  ensemble_trace, 
  data.frame(time = TWICE,count = ensemble_trace[nrow(ensemble_trace),"count"], event = "end")
)

# structure for Plotting
ensemble_trace_df <- data.frame(x0=-1,y0=-1,x1=-1,y1=-1,id=-1)
ensemble_trace_vert_df <- data.frame(x0=-1,y0=-1,x1=-1,y1=-1,id=-1)
for(i in 1:(nrow(ensemble_trace)-1)){
  
  ensemble_trace_df <- rbind(
    ensemble_trace_df,
    data.frame(
      x0 = ensemble_trace[i,"time"],
      x1 = ensemble_trace[i+1,"time"],
      y0 = ensemble_trace[i,"count"],
      y1 = ensemble_trace[i,"count"],
      id = i
    )
  )
  
  if(i < (nrow(ensemble_trace)-1)){
    ensemble_trace_vert_df <- rbind(
      ensemble_trace_vert_df,
      data.frame(
        x0 = ensemble_trace[i+1,"time"],
        x1 = ensemble_trace[i+1,"time"],
        y0 = ensemble_trace[i,"count"],
        y1 = ensemble_trace[i+1,"count"],
        id = i
      )
    )  
  }
  
}
ensemble_trace_df <- ensemble_trace_df[-1,]
ensemble_trace_vert_df <- ensemble_trace_vert_df[-1,]

# plots
plot_indiv <- ggplot(data = indiv_traces) +
  geom_segment(data = indiv_traces,aes(y=id,x=enter,yend=id,xend=leave-0.05),size=1.025) +
  geom_point(aes(x=enter,y=id),shape=16,size=4)+
  geom_point(aes(x=leave,y=id),shape=1,size=4)+
  scale_y_continuous(breaks=indiv_traces$id,labels=as.character(indiv_traces$id)) +
  labs(x = "Time",title = "Individual Life Histories") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),axis.text.y = element_text(size=14),axis.title.y = element_blank(),
    axis.title.x = element_text(size=16),title = element_text(size=18)
  )

ggsave(filename = here::here("figs/indiv_hist.tiff"),plot = plot_indiv,device = "tiff",height = 6,width = 8,compression = "lzw")

plot_ensemble <- ggplot(data = ensemble_trace_df) + 
  geom_segment(aes(y=y0,x=x0,yend=y1,xend=x1,group=id),size=1.025) +
  geom_segment(data=ensemble_trace_vert_df,aes(x=x0,xend=x1,y=y0,yend=y1,group=id),linetype=2,alpha=0.5) +
  geom_point(aes(x=x0,y=y0),shape=16,size=4) +
  geom_point(aes(x=x1,y=y1),shape=1,size=4) +
  geom_text(data=ensemble_trace[!ensemble_trace$event %in% c("begin","end"),], aes(x=time,y=count-0.05,label=event),vjust = 1,hjust = 1,size=2) +
  labs(x = "Time",title = "Ensemble Trajectory") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),axis.text.y = element_text(size=14),axis.title.y = element_blank(),
    axis.title.x = element_text(size=16),title = element_text(size=18)
  )

ggsave(filename = here::here("figs/ensemble_hist.tiff"),plot = plot_ensemble,device = "tiff",height = 6,width = 8,compression = "lzw")
