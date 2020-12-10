




make_mosquito <- function(tnow, parameters){
  m <- list(
    thist = tnow,
    shist = "S",
    tdie = tnow + rexp(n = 1,rate = parameters[["g"]]),
    H2M_bite = NaN
  )
  m$tnext <- m$tdie
}

