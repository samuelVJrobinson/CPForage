#Create and optimize currency in a world where TRANSFER foragers has been added to nest(s) WHICHNEST
makeBest=function(scenario,whichNest=NA,parallel=F,cluster=NA){

  if(parallel&&is.na(cluster)) stop('Cluster not specified')
  #Not used currently, but could be used to restrict which nest to add foragers to (currently added to all). Essentially represent "worst case scenario".
  if(is.na(whichNest)) stop('Nest # not specified')

  for(i in 1:length(scenario$nests)){ #Add foragers to each cell in each nest
    transfer=with(scenario$nests[[i]],steps[stepNum]) #Number of foragers to add
    scenario$nests[[i]]$n=scenario$nests[[i]]$n+transfer #Adds TRANSFER foragers to every cell
  }
  use=matrix(T,nrow=nrow(scenario$nests[[i]]$n),ncol=ncol(scenario$nests[[i]]$n)) #All cells
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
