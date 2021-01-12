
# NH = size of human pop
# X = prevalence in humans
# g = 1/avg lifespan
# a = human bloodfeeding rate
# b = mosy -> human transmission efficiency
# c = human -> mosy transmission efficiency
# r = 1/avg duration of infection in humans
# EIP = length of EIP
# LEP = duration in liver
calc_equilibrium <- function(NH,X,g=1/12,a=(0.9 * 1/3),b=0.55,c=0.15,r=1/200,EIP=10,LEP=7){

  P <- exp(-g*EIP)
  SH <- NH*(1-X)
  IH <- NH-SH

  # equilibrium IV
  IV <- (r*IH*(SH + IH)) / (a * b * SH)

  # prob mosquito becomes infected and survives to become infectious
  p_Minf <- (a*c*X) / ((a*c*X) + g)
  p_Minf <- p_Minf * P

  # total num mosy
  NV <- IV / p_Minf

  # solve for lambda
  lambda <- NV*g

  # calc EV
  EV <- (1/g) * ( ((g*IV)/P)   - (g*IV) )

  # susceptible vectors
  SV <- NV - IV - EV

  list(
    "y0" = c(SH=SH,IH=IH,SV=SV,EV=EV,IV=IV),
    "parameters" = c(g=g,b=b,a=a,c=c,r=r,EIP=EIP,LEP=LEP,SV0=SV,IV0=IV,SH0=SH,IH0=IH,lambda=lambda)
  )
}

# # for deSolve::dede
# derivs <- function(t,y,pars) {
#
#   with(as.list(pars),{
#     SH <- y[1]
#     IH <- y[2]
#     SV <- y[3]
#     IV <- y[4]
#
#     if(t <= EIP){
#       SV_eip <- SV0
#
#       SH_eip <- SH0
#       IH_eip <- IH0
#     } else {
#       SV_eip <- lagvalue(t - EIP,3)
#
#       SH_eip <- lagvalue(t - EIP,1)
#       IH_eip <- lagvalue(t - EIP,2)
#     }
#
#     if(t <= LEP){
#       SH_lep <- SH0
#       IH_lep <- IH0
#
#       IV_lep <- IV0
#     } else {
#       SH_lep <- lagvalue(t - LEP,1)
#       IH_lep <- lagvalue(t - LEP,2)
#
#       IV_lep <- lagvalue(t - LEP,4)
#     }
#
#     X_lep <- IH_lep / (SH_lep + IH_lep)
#     X_eip <- IH_eip / (SH_eip + IH_eip)
#     X <- IH / (SH + IH)
#
#     dSH <- (r * IH) - (a * b * (1 - X) * IV)
#     dIH <- (a * b * IV_lep * (1 - X_lep)) - (r * IH)
#
#     dSV <- lambda - (a * c * X * SV) - (g * SV)
#     dIV <- (a * c * X_eip * SV_eip * exp(-g * EIP)) - (g * IV)
#
#     return(list(c(dSH,dIH,dSV,dIV)))
#   })
# }

derivs_odin <- odin::odin({

  # lags
  SH_lep <- delay(SH,LEP)
  IH_lep <- delay(IH,LEP)
  IV_lep <- delay(IV,LEP)

  SH_eip <- delay(SH,EIP)
  IH_eip <- delay(IH,EIP)
  SV_eip <- delay(SV,EIP)

  # prevalence in human population
  X <- IH / (SH + IH)
  X_lep <- IH_lep / (SH_lep + IH_lep)
  X_eip <- IH_eip / (SH_eip + IH_eip)

  deriv(SH) <- (r * IH) - (a * b * (1 - X) * IV)
  deriv(IH) <- (a * b * IV_lep * (1 - X_lep)) - (r * IH)

  deriv(SV) <- lambda - (a * c * X * SV) - (g * SV)
  deriv(IV) <- (a * c * X_eip * SV_eip * exp(-g * EIP)) - (g * IV)

  initial(SH) <- SH0
  initial(IH) <- IH0
  initial(SV) <- SV0
  initial(IV) <- IV0

  # parameters and initial conditions
  SH0 <- user()
  IH0 <- user()
  SV0 <- user()
  IV0 <- user()

  g <- user()
  b <- user()
  a <- user()
  c <- user()
  r <- user()
  EIP <- user()
  LEP <- user()
  lambda <- user()
})
