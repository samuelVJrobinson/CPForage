#Create and optimize currency in a world where TRANSFER foragers has been added to nest(s) WHICHNEST
makeBest=function(nests,world,whichNest=NA,parallel=F,cluster=NA){
  if(parallel&&is.na(cluster)) stop('Cluster not specified')
  #Not used currently, but could be used to restrict which nest to add foragers to (currently added to all). Essentially represent "worst case scenario".
  if(is.na(whichNest)) stop('Nest # not specified')
  bestNests=nests
  bestWorld=world
  for(i in 1:length(bestNests)){ #Add foragers to each cell in each nest
    transfer=bestNests[[i]]$steps[bestNests[[i]]$stepNum] #Number of foragers to add
    bestNests[[i]]$n=bestNests[[i]]$n+transfer #Adds TRANSFER foragers to every cell
  }
  use=matrix(T,nrow=nrow(bestNests[[1]]$n),ncol=ncol(bestNests[[1]]$n)) #All cells
  #Calculates S,L, and Curr values for every nest
  if(parallel){
    temp=parLapply(cluster,which(use),optimLoadCurr,nests=bestNests,world=bestWorld)
  } else {
    temp=lapply(which(use),optimLoadCurr,nests=bestNests,world=bestWorld)
  }
  bestWorld$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-value
  for(u in 1:length(temp)){ #For each cell processed
    for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
      bestNests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
      bestNests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns curr
    }
  }
  return(list('bestNests'=bestNests,'bestWorld'=bestWorld))
}
