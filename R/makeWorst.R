#'World with -TRANSFER less foragers
#'
#Create and optimize currency in a world where -TRANSFER foragers has been taken from nest(s) WHICHNEST
#'
#'@param scenario Nest structure and world structure (scenario)
#'@param whichNest Which nest should have TRANSFER foragers taken away?
#'@param parallel Should computation be done in parallel?
#'@param cluster If parallel, which cluster should be used

#'@return List of nests and world structure (scenario)
#'
#'@examples

makeWorst=function(scenario,whichNest=NA,parallel=F,cluster=NA) {
  if(parallel&&is.na(cluster)) stop('Cluster not specified')
  #Not used currently, but could be used to restrict which nest to add foragers to (currently added to all)
  if(is.na(whichNest)) stop('Nest # not specified')

  for(i in whichNest){ #Calculates load size and currency for each nest
    transfer=scenario$nests[[i]]$steps[scenario$nests[[i]]$stepNum] #Number of foragers to remove
    use=scenario$nests[[i]]$n>=transfer #Occupied cells with at least TRANSFER foragers in them
    scenario$nests[[i]]$n[use]=scenario$nests[[i]]$n[use]-transfer #Subtracts TRANSFER forager from every occupied cell
  }
  #Calculates S,L, and Curr values for every nest
  if(parallel){
    temp=parLapply(cluster,which(use),optimLoadCurr,scenario=scenario)
  } else {
    temp=lapply(which(use),optimLoadCurr,scenario=scenario)
  }
  scenario$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-value
  for(u in 1:length(temp)){ #For each cell processed
    for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
      scenario$nests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
      scenario$nests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns curr
    }
  }
  return(scenario)
}
