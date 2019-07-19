#'@title Create scenario with +TRANSFER more foragers
#'
#'@description Create and optimize currency in a world where -\code{transfer} foragers has been taken from nest(s) \code{whichnest}. Called within \code{forageMod}
#'
#'@param scenario Nest structure and world structure (scenario)
#'@param whichNest Which nest should have TRANSFER foragers taken away?
#'@param parallel Should computation be done in parallel?
#'@param cluster If parallel, which cluster should be used

#'@return List of nests and world structure (scenario)
#'
#'@examples
#'makeBest(base,whichNest=1)
makeBest=function(scenario,whichNest=NA,parallel=FALSE,cluster=NA){
  if(parallel&&is.na(cluster)) stop('Cluster not specified')
  #Not used currently, but could be used to restrict which nest to add foragers to (currently added to all). Essentially represent "worst case scenario", where everyone moves to a given cell.
  if(is.na(whichNest)) stop('Nest # not specified')

  for(i in 1:length(scenario$nests)){ #Add foragers to each cell in each nest
    transfer=with(scenario$nests[[i]],steps[stepNum]) #Number of foragers to add
    scenario$nests[[i]]$n=scenario$nests[[i]]$n+transfer #Adds TRANSFER foragers to every cell
  }
  use=matrix(TRUE,nrow=nrow(scenario$nests[[i]]$n),ncol=ncol(scenario$nests[[i]]$n)) #All cells
  #Calculates S,L, and Curr values for every nest
  if(parallel){
    temp=parLapply(cluster,which(use),optimLoadCurr,scenario=scenario)
  } else {
    temp=lapply(which(use),optimLoadCurr,scenario=scenario)
  }
  scenario$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-value
  scenario$nests[[1]]$L[which(use)]=sapply(temp,function(x) x$optimL) #Assigns L
  scenario$nests[[1]]$curr[which(use)]=sapply(temp,function(x) x$optimCurr) #Assigns curr

  return(scenario)
}
