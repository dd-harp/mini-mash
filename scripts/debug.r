rm(list=ls());gc()
dev.off()

library(here)
library(data.table)
library(ggplot2)
library(deSolve)
library(Rcpp)
library(stocheulerABM)

# source(here::here("r-src/deterministic-exposed.R"))
# 
# NH <- 1e3
# X <- 0.3
# 
# IC <- calc_equilibrium(NH = NH,X = X)
# 
# 
# 
# tmax <- 50
# 
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
#   EH0 = 250,
#   EV0 = 0,
#   SH0 = IC$y0[["SH"]],
#   SV0 = IC$y0[["SV"]]
# )
# 
# dde_out <- dde_odin$run(t = 0:tmax)
# 
# dde_out <- as.data.table(dde_out)
# data.table::setnames(x = dde_out,old = "t",new = "time")
# dde_out <- data.table::melt(dde_out,id.vars = "time")
# 
# ggplot(data = dde_out) +
#   geom_line(aes(x=time,y=value,color=variable),size=1.05) +
#   theme_bw()
# 
# dde_out[variable=="EH",]
# dde_out[variable=="EV",]


# the blood meal example
tmin <- 0
tmax <- 5

set.seed(342)

St <- c(tmin,sort(runif(n = 2,min = 0,max = 5)),tmax)
Sv <- c(100,rpois(n = 2,lambda = 100),100)

Xt <- c(tmin,sort(runif(n = 1,min = 0,max = 5)),tmax)
Xv <- c(0.5,rbeta(n = 1,shape1 = 10,shape2 = 10),0.5)



a <- 0.9 * 1.3
c <- 0.15

LambdaAnaly <- S_traj[[2]][1] * 100 * a * c * 0.5
LambdaAnaly <- LambdaAnaly +( (4.8778720 -  1.103043 ) * a * c * 90 * 0.5 )
LambdaAnaly <- LambdaAnaly +( ( 4.958862 - 4.8778720 ) * a * c * 90 * 0.5649469 )
LambdaAnaly <- LambdaAnaly +( ( 5 - 4.958862 ) * a * c * 112 * 0.5649469 )

# algorithm
S_traj <- mapply(FUN = function(x,y){
  c(x,y)
},x=St,y=Sv,SIMPLIFY = F)

X_traj <- mapply(FUN = function(x,y){
  c(x,y)
},x=Xt,y=Xv,SIMPLIFY = F)

lambda <- 0
t0 <- tmin

t0_S <- S_traj[[1]]
S_traj <- S_traj[-1]

t0_X <- X_traj[[1]]
X_traj <- X_traj[-1]

t1_S <- S_traj[[1]]
S_traj <- S_traj[-1]

t1_X <- X_traj[[1]]
X_traj <- X_traj[-1]

t1_times <- c(S=t1_S[1],X=t1_X[1])
i <- 0
ddt <- 0
# while( !( length(S_traj)==0 & length(X_traj)==0 )) {
while( all(is.finite(t1_times)) )   {
  
  mu <- which.min(t1_times)
  t1 <- t1_times[mu]
  dt <- t1 - t0
  ddt <- ddt + dt
  
  lambda <- lambda + (a * c * t0_X[2] * t0_S[2] * dt)
    
  t0 <- t1
  if(mu==1){
    t0_S <- t1_S
    if( length(S_traj)!=0 ){
      t1_S <- S_traj[[1]]
      S_traj <- S_traj[-1]
      t1_times[1] <- t1_S[1]
    } else {
      t1_times[1] <- Inf
    }
  } else {
    t0_X <- t1_X
    if( length(X_traj)!=0 ){
      t1_X <- X_traj[[1]]
      X_traj <- X_traj[-1]
      t1_times[2] <- t1_X[1]
    } else {
      t1_times[2] <- Inf
    }
  }
  i <- i +1
}





