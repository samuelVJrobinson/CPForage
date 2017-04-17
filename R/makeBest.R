#Create and optimize currency in a world where TRANSFER foragers has been added to nest(s) WHICHNEST
makeBest=function(nests,world,whichNest=NA,parallel=F,cluster=NA){
  if(parallel&&is.na(cluster)) stop('Cluster not specified')
  if(is.na(whichNest)) stop('Nest # not specified')
  bestNests=nests
  bestWorld=world
  for(i in 1:length(bestNests)){ #Calculates initial loading rate, load size, and currency for each nest
    transfer=bestNests[[i]]$steps[bestNests[[i]]$stepNum] #Number of foragers to remove
    bestNests[[i]]$n=bestNests[[i]]$n+transfer #Adds TRANSFER foragers to every cell
    use=matrix(T,nrow=nrow(bestNests[[i]]$n),ncol=ncol(bestNests[[i]]$n)) #All cells
    #Calculates S,L, and Curr values for every nest
    if(parallel){
      temp=parLapply(c1=cluster,which(use),optimLoadCurr,nests=bestNests,world=bestWorld)
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
  }
  return(list('bestNests'=bestNests,'bestWorld'=bestWorld))
}
