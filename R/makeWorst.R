#Create and optimize currency in a world where -TRANSFER foragers has been taken from nest(s) WHICHNEST
makeWorst=function(nests,world,whichNest=NA,parallel=F,cluster=NA) {
  if(parallel&&is.na(cluster)) stop('Cluster not specified')
  if(is.na(whichNest)) stop('Nest # not specified')
  worstNests=nests
  worstWorld=world
  for(i in whichNest){ #Calculates load size and currency for each nest
    transfer=worstNests[[i]]$steps[worstNests[[i]]$stepNum] #Number of foragers to remove
    use=worstNests[[i]]$n>=transfer #Occupied cells with at least TRANSFER foragers in them
    worstNests[[i]]$n[use]=worstNests[[i]]$n[use]-transfer #Subtracts TRANSFER forager from every occupied cell
  }
  occupied=do.call('+',lapply(nests, function(x) x$n>0))>0 #Cells with >0 foragers
  #Calculates S,L, and Curr values for every nest
  if(parallel){
    temp=parLapply(cluster,which(use),optimLoadCurr,nests=worstNests,world=worstWorld)
  } else {
    temp=lapply(which(use),optimLoadCurr,nests=worstNests,world=worstWorld)
  }
  worstWorld$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-value
  for(u in 1:length(temp)){ #For each cell processed
    for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
      worstNests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
      worstNests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns curr
    }
  }
  #Sets values in unoccupied cells to NA (nothing can be moved from them )
  occupied=worstNests[[i]]$n>0
  worstNests[[i]]$L[!occupied]=worstNests[[i]]$curr[!occupied]=NA
  return(list(worstNests=worstNests,worstWorld=worstWorld))
}
