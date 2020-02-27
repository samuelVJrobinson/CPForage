#'@title Create scenario with -TRANSFER less foragers
#'
#'@description Create and optimize currency in a world where -TRANSFER foragers has been taken from nest(s) WHICHNEST. Called within \code{forageMod}
#'
#'@param scenario Nest structure and world structure (scenario)
#'@param parallel Should computation be done in parallel?
#'@param cluster If parallel, which cluster should be used

#'@return List of nests and world structure (scenario)
#'
#'@examples
#'makeWorst(base,parallel=T,cluster=cluster)

makeWorst=function(scenario,parallel=F,cluster=NA) {
  if(parallel&&is.na(cluster)) stop('Cluster not specified')

  transfer <- scenario$nests$steps[scenario$nests$stepNum] #Number of foragers to remove
  use <- scenario$nests$n>=transfer #Occupied cells with at least TRANSFER foragers in them
  scenario$nests$n[use] <- scenario$nests$n[use] - transfer #Subtracts TRANSFER forager from every occupied cell

  #Calculates S,L, and Curr values for every nest
  if(parallel){
    temp <- parLapply(cluster,which(use),optimLoadCurr,scenario=scenario)
  } else {
    temp <- lapply(which(use),optimLoadCurr,scenario=scenario)
  }
  scenario$world$S[which(use)] <- sapply(temp,function(x) x$S) #Assigns S-value
  scenario$nests$L[which(use)] <- sapply(temp,function(x) x$optimL) #Assigns L
  scenario$nests$curr[which(use)] <- sapply(temp,function(x) x$optimCurr) #Assigns curr

  return(scenario)
}
