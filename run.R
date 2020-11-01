library(here)
library(data.table)
library(ggplot2)
library(Rcpp)
library(stocheulerABM)

Rcpp::sourceCpp(here::here("r-src/aggregated.cpp"))

params <- c("lambda" = 5, "g" = 1/12, "b" = 0.55, "c" = 0.15, "a" = 0.9 * 1/3, "r" = 1/200, "EIP" = 5, "LEP" = 7)
tmax <- 365*50

# out <- pfsim_aggregated(
#   tmax = tmax,
#   SH = 500,
#   IH = 10,
#   SV = 100,
#   IV = 5,
#   parameters = params,
#   verbose = T
# )
# 
# out <- data.table::as.data.table(discretise(out = out,dt = 1))
# out <- data.table::melt(out,id.vars = "time",measure.vars = c("SH","IH","SV","IV"))
# 
# ggplot(data = out) +
#   geom_line(aes(x=time,y=value,color=variable)) +
#   theme_bw()





SH0 <- 100
IH0 <- 50


calc_equilibrium <- function(parameters, SH, IH){
  with(as.list(parameters),{
    X <- (IH / (IH + SH))
    IV <- (r * IH) / (a * b * (1 - X))
    SV <- (g * IV) / (a * c * X * exp(-g*EIP))
    lambda <- SV * ((a * c * X) + g)
    return(c("SV"=SV,"IV"=IV,"lambda"=lambda))
  })
}

eq <- calc_equilibrium(parameters = params,SH = 500,IH = 100)
pars1 <- params
pars1["lambda"] <- eq["lambda"]

out <- pfsim_aggregated(
  tmax = tmax,
  SH = SH0,
  IH = IH0,
  SV = floor(eq["SV"]),
  IV = floor(eq["IV"]),
  parameters = pars1,
  verbose = T
)

out <- data.table::as.data.table(discretise(out = out,dt = 1))
out <- data.table::melt(out,id.vars = "time",measure.vars = c("SH","IH","SV","IV"))

ggplot(data = out) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()

library(deSolve)


derivs <- function(t,y,pars) {
  
  with(as.list(pars),{
    SH <- y[1]
    IH <- y[2]
    SV <- y[3]
    IV <- y[4]
    
    if(t <= EIP){
      SV_eip <- SV0
      
      SH_eip <- SH0
      IH_eip <- IH0
    } else {
      SV_eip <- lagvalue(t - EIP,3)
      
      SH_eip <- lagvalue(t - EIP,1)
      IH_eip <- lagvalue(t - EIP,2)
    }
    
    if(t <= LEP){
      SH_lep <- SH0
      IH_lep <- IH0
      
      IV_lep <- IV0
    } else {
      SH_lep <- lagvalue(t - LEP,1)
      IH_lep <- lagvalue(t - LEP,2)
      
      IV_lep <- lagvalue(t - LEP,4)
    }
    
    X_lep <- IH_lep / (SH_lep + IH_lep)
    X_eip <- IH_eip / (SH_eip + IH_eip)
    X <- IH / (SH + IH)
    
    dSH <- (r * IH) - (a * b * (1 - X) * IV)
    dIH <- (a * b * IV_lep * (1 - X_lep)) - (r * IH)
    
    dSV <- lambda - (a * c * X * SV) - (g * SV)
    dIV <- (a * c * X_eip * SV_eip * exp(-g * EIP)) - (g * IV)
    
    return(list(c(dSH,dIH,dSV,dIV)))
  })
}


desolve_pars <- pars1
desolve_pars[["SH0"]] <- SH0
desolve_pars[["IH0"]] <- IH0
desolve_pars[["SV0"]] <- eq[["SV"]]
desolve_pars[["IV0"]] <- eq[["IV"]]

y0 <- c("SH"=SH0,"IH"=IH0,"SV"=eq[["SV"]],"IV"=eq[["IV"]])

dde_out <- dede(y = y0,times = 0:tmax,func = derivs,parms = desolve_pars)
plot(dde_out)

dde_out <- as.data.table(dde_out)
dde_out <- data.table::melt(dde_out,id.vars = "time",measure.vars = c("SH","IH","SV","IV"))

out[, ("type") := "stochastic"]
dde_out[, ("type") := "dde"]

ggplot(data = rbind(out,dde_out)) +
  geom_line(aes(x=time,y=value,color=variable,linetype=type),alpha=0.85) +
  theme_bw()



IHf <- dde_out[variable=="IH" & time == max(time),"value"][[1]]
SHf <- dde_out[variable=="SH" & time == 18250,"value"][[1]]
SVf <- dde_out[variable=="SV" & time == 18250,"value"][[1]]
IVf <- dde_out[variable=="IV" & time == 18250,"value"][[1]]
Xf <- (IHf / (IHf + SHf))

# ok for IV

(desolve_pars[["r"]]*IH0) 

(desolve_pars[["r"]]*IHf) / (desolve_pars[["a"]] * desolve_pars[["b"]]  * (1 - Xf))

# SV
(desolve_pars[["g"]] * IVf) / (desolve_pars[["a"]] * desolve_pars[["c"]] * Xf * exp(-desolve_pars[["g"]] * desolve_pars[["EIP"]]))

# lambda
SVf * ((desolve_pars[["a"]] * desolve_pars[["c"]] * Xf) + desolve_pars[["g"]])






# try again
g <- params[["g"]]
b <- params[["b"]]
c <- params[["c"]]
a <- params[["a"]]
r <- params[["r"]]
EIP <- params[["EIP"]]
LEP <- params[["LEP"]]
P <- exp(-g*EIP)

(r*IH0) / (a * b * (SH0 / (SH0 + IH0)))
