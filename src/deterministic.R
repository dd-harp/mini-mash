
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

  # P(survive EIP)
  P <- exp(-g*EIP)

  # infectious human population
  IH <- NH*X

  # incubating human population
  EH <- IH*LEP*r

  # susceptible human population
  SH <- NH - IH - EH

  # equilibrium IV
  IV <-  -(IH*NH*r) / (a*b * (IH - NH + (IH*LEP*r)))

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
    "y0" = c(SH=SH,EH=EH,IH=IH,SV=SV,EV=EV,IV=IV),
    "parameters" = c(g=g,b=b,a=a,c=c,r=r,EIP=EIP,LEP=LEP,SV0=SV,IV0=IV,SH0=SH,IH0=IH,lambda=lambda)
  )
}


derivs_odin <- odin::odin({

  # lags
  SH_lep <- delay(SH,LEP)
  EH_lep <- delay(EH,LEP)
  IH_lep <- delay(IH,LEP)
  IV_lep <- delay(IV,LEP)

  SH_eip <- delay(SH,EIP)
  EH_eip <- delay(EH,EIP)
  IH_eip <- delay(IH,EIP)
  SV_eip <- delay(SV,EIP)

  Z <- SH / (SH + EH + IH)
  Z_lep <- SH_lep / (SH_lep + EH_lep + IH_lep)

  X <- IH / (SH + EH + IH)
  X_eip <- IH_eip / (SH_eip + EH_eip + IH_eip)

  deriv(SH) <- (r * IH) - (a * b * Z * IV)
  deriv(EH) <- (a * b * Z * IV) - (a * b * IV_lep * Z_lep)
  deriv(IH) <- (a * b * IV_lep * Z_lep) - (r * IH)

  deriv(SV) <- lambda - (a * c * X * SV) - (g * SV)
  deriv(EV) <- (a * c * X * SV) - (a * c * X_eip * SV_eip * exp(-g * EIP)) - (g * EV)
  deriv(IV) <- (a * c * X_eip * SV_eip * exp(-g * EIP)) - (g * IV)

  initial(SH) <- SH0
  initial(EH) <- EH0
  initial(IH) <- IH0
  initial(SV) <- SV0
  initial(EV) <- EV0
  initial(IV) <- IV0

  # parameters and initial conditions
  SH0 <- user()
  EH0 <- user()
  IH0 <- user()
  SV0 <- user()
  EV0 <- user()
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
