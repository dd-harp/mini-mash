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


#  --------------------------------------------------------------------------------
# test abm
# --------------------------------------------------------------------------------

# Rcpp::sourceCpp(here::here("r-src/exposed/disaggregated-abm.cpp"),showOutput = FALSE)
# 
# abm_IC <- round(IC$y0)
# SH <- abm_IC[["SH"]]
# IH <- sum(abm_IC[c("EH","IH")])
# SV <- abm_IC[["SV"]]
# IV <- sum(abm_IC[c("EV","IV")])
# 
# set.seed(342342L)
# abmout <- run_abm(SV = SV,IV = IV,SH = SH,IH = IH,parameters = IC$parameters,dt = 5,tmax = 5e3)
# hout <- abmout$human
# mout <- abmout$mosy
# 
# init_ix <- which(hout[,1]<=0)
# init <- colSums(hout[init_ix,-1])
# hout <- hout[-init_ix,]
# 
# hhist <- matrix(data = 0,nrow = nrow(hout)+1,ncol = 4)
# hhist[1,2:4] <- init
# hhist[2:nrow(hhist),1] <- hout[,1]
# for(i in 1:nrow(hout)){
#   hhist[i+1,2:4] <- hhist[i,2:4] + hout[i,2:4]
# }
# 
# init_ix <- which(mout[,1] <= 0)
# init <- colSums(mout[init_ix,-1])
# mout <- mout[-init_ix,]
# 
# mhist <- matrix(data = 0,nrow = nrow(mout)+1,ncol = 4)
# mhist[1,2:4] <- init
# mhist[2:nrow(mhist),1] <- mout[,1]
# for(i in 1:nrow(mout)){
#   mhist[i+1,2:4] <- mhist[i,2:4] + mout[i,2:4]
# }

# # compare
# colMeans(hhist[,2:4])
# unname(IC$y0[1:3])
# 
# colMeans(mhist[,2:4])
# unname(IC$y0[c("SV","EV","IV")])

abm_IC <- round(IC$y0)
SH <- abm_IC[["SH"]]
IH <- sum(abm_IC[c("EH","IH")])
SV <- abm_IC[["SV"]]
IV <- sum(abm_IC[c("EV","IV")])

Rcpp::sourceCpp(here::here("r-src/exposed/disaggregated-abm.cpp"),showOutput = FALSE)

# set.seed(342342L)
abmoutsum <- run_abm_sumout(SV = SV,IV = IV,SH = SH,IH = IH,parameters = IC$parameters,dt = 5,tmax = 1e4)

colMeans(abmoutsum$human[,2:4])
IC$y0[1:3]

colMeans(abmoutsum$mosy[,2:4])
IC$y0[c("SV","EV","IV")]

#  --------------------------------------------------------------------------------
# test humans
# --------------------------------------------------------------------------------

Rcpp::sourceCpp(here::here("r-src/exposed/test-human-abm.cpp"),showOutput = FALSE)


hum_IC <- round(IC$y0[1:3] )

hout <- test_humans(SH = hum_IC[1],IH = sum(hum_IC[2:3]),IV = IC$y0[["IV"]],parameters = IC$parameters,dt = 5,tmax = 5e3)

# state_df <- as.data.frame(hout)
# inf_per <- NULL
# exp_per <- NULL
# for(i in 1:(nrow(state_df)-1)){
#   stateseq <- state_df[i:(i+1),"states"]
#   if(all(stateseq == c("E","I"))){
#     exp_per <- c(exp_per, (state_df[(i+1),"times"] - state_df[i,"times"]) )
#   }
#   if(all(stateseq == c("I","S"))){
#     inf_per <- c(inf_per, (state_df[(i+1),"times"] - state_df[i,"times"]) )
#   }
# }
# 
# hist(inf_per[inf_per>0])
# mean(inf_per[inf_per>0])
# 
# state_dt <- as.data.table(state_df)
sout <- hout$agghist

init_ix <- which(sout[,1]<=0)
init <- colSums(sout[init_ix,-1])
sout <- sout[-init_ix,]

shist <- matrix(data = 0,nrow = nrow(sout)+1,ncol = 4)
shist[1,2:4] <- init
shist[2:nrow(shist),1] <- sout[,1]
for(i in 1:nrow(sout)){
  shist[i+1,2:4] <- shist[i,2:4] + sout[i,2:4]
}
matplot(x = shist[,1],y = shist[,2:4],lty=1,type="l")
colMeans(shist[,2:4])
IC$y0[1:3]


#  --------------------------------------------------------------------------------
# test mosy
# --------------------------------------------------------------------------------

Rcpp::sourceCpp(here::here("r-src/exposed/test-mosy-abm.cpp"),showOutput = FALSE)

mosy_IC <- round(IC$y0[c("SV","EV","IV")])
SV <- mosy_IC[1]
IV <- sum(mosy_IC[2:3])

IH <- IC$y0[["IH"]]
NH <- sum(IC$y0[c("SH","EH","IH")])

mout <- test_mosy(
  SV = SV,IV = IV,IH = IH,NH = NH,
  parameters = IC$parameters,dt = 5,tmax = 5e3
)

init_ix <- which(mout[,1] <= 0)
init <- colSums(mout[init_ix,-1])
mout <- mout[-init_ix,]

mhist <- matrix(data = 0,nrow = nrow(mout)+1,ncol = 4)
mhist[1,2:4] <- init
mhist[2:nrow(mhist),1] <- mout[,1]
for(i in 1:nrow(mout)){
  mhist[i+1,2:4] <- mhist[i,2:4] + mout[i,2:4]
}

matplot(x = mhist[,1],y = mhist[,2:4],type="l",lty=1)

colMeans(mhist[,2:4])
IC$y0[c("SV","EV","IV")]
