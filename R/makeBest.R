#'@title Create scenario with +TRANSFER more foragers
#'
#'@description Create and optimize currency in a world where -\code{transfer} foragers has been taken from nest(s) \code{whichnest}. Called within \code{forageMod}
#'
#'@param scenario Nest structure and world structure (scenario)
#'@param parallel Should computation be done in parallel?
#'@param cluster If parallel, which cluster should be used

#'@return List of nests and world structure (scenario)
#'
#'@examples
#'makeBest(base,parallel=TRUE,cluster=cluster)
makeBest <- function(scenario,parallel=FALSE,cluster=NA){
  if(parallel&&is.na(cluster)) stop('Cluster not specified')

  transfer <- with(scenario$nests,steps[stepNum]) #Number of foragers to add
  scenario$nests$n <- scenario$nests$n + transfer #Adds TRANSFER foragers to every cell

  use <- matrix(T,nrow=nrow(scenario$nests$n),ncol=ncol(scenario$nests$n)) #All cells
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
