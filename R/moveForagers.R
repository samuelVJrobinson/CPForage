moveForagers=function(scenarioSet,i,moves){ #Function to move foragers and update scenarios in nest i, given scenarioSet and moves list
  if(!moves$move){
    warning('Move not required')
    return(scenarioSet)
  }
  #Unpack moves list
  worst=moves$from
  best=moves$to

  #Unpack scenarioSet - currently works for single-forager case, or multi-forager situation ignoring other nests
  bestScen=scenarioSet$best
  baseScen=scenarioSet$base
  sol=baseScen$nests[[i]]$sol #Solitary foraging?
  if(!sol) worstScen=scenarioSet$worst
  transfer=baseScen$nests[[i]]$steps[baseScen$nests[[i]]$stepNum] #Transfer number

  if(sol){ #If foragers are solitary

    #Step 1: Copy values from bestNests,bestWorld to nests,world

    #Overwrite TRANSFER foragers to the best cell
    baseScen$nests[[i]]$n[best]=bestScen$nests[[i]]$n[best]
    #Overwrite optimal Load and Currency from the best cell
    #NOTE: should be done for all nests, since more than nest i may have changed
    baseScen$nests[[i]]$L[best]=bestScen$nests[[i]]$L[best] #Load size
    baseScen$nests[[i]]$curr[best]=bestScen$nests[[i]]$curr[best] #Currency
    baseScen$world$S[best]=bestScen$world$S[best] #S-value

    #Step 2: Subtract foragers from nests, and bestNests

    #Move TRANSFER foragers out of worst cell in NESTS
    baseScen$nests[[i]]$n[worst]=baseScen$nests[[i]]$n[worst]-transfer
    #Move TRANSFER foragers out of (and into) worst (and best) cells in BESTNESTS
    bestScen$nests[[i]]$n[best]=bestScen$nests[[i]]$n[best]+transfer #Adds new foragers from best cell
    bestScen$nests[[i]]$n[worst]=bestScen$nests[[i]]$n[worst]-transfer #Subtracts foragers to worst cell

    #Step 3: Update values for nests and bestNests

    #For nests:
    temp=lapply(which(worst),optimLoadCurr,scenario=baseScen)
    baseScen$world$S[which(worst)]=sapply(temp,function(x) x$S) #Assigns S-values
    for(u in 1:length(temp)){ #For each cell processed
      for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
        baseScen$nests[[name]][['L']][which(worst)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
        baseScen$nests[[name]][['curr']][which(worst)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns currency
      }
    }

    #For bestNests:
    use=worst|best #Location of worst and best cells
    temp=lapply(which(use),optimLoadCurr,scenario=bestScen)
    bestScen$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
    for(u in 1:length(temp)){ #For each cell processed
      for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
        bestScen$nests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
        bestScen$nests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns currency
      }
    }

  } else { #If foragers are social

    #Step 1a: Copy values from bestNests,bestWorld to nests,world

    #Copy TRANSFER foragers to best cell
    baseScen$nests[[i]]$n[best]=bestScen$nests[[i]]$n[best] #Overwrite forager number for best cell
    baseScen$nests[[i]]$L[best]=bestScen$nests[[i]]$L[best] #Load size
    baseScen$nests[[i]]$curr[best]=bestScen$nests[[i]]$curr[best] #Currency
    baseScen$world$S[best]=bestScen$world$S[best] #S-value

    #Step 1b: Copy values from worstNests,worstWorld to nests,world

    #Copy TRANSFER foragers from worst cell
    baseScen$nests[[i]]$n[worst]=worstScen$nests[[i]]$n[worst] #Overwrite forager number for worst cell
    baseScen$nests[[i]]$L[worst]=worstScen$nests[[i]]$L[worst] #Load size
    baseScen$nests[[i]]$curr[worst]=worstScen$nests[[i]]$curr[worst] #Currency
    baseScen$world$S[worst]=worstScen$world$S[worst] #S-value


    #Step 2: Subtract foragers from worstNests, and bestNests

    #Move TRANSFER foragers out of (and into) worst (and best) cells in WORSTNESTS
    worstScen$nests[[i]]$n[best]=worstScen$nests[[i]]$n[best]+transfer #Adds new foragers to best cell
    #Subtracts foragers from worst cell - if <= transfer, sets foragers to 0 (won't be considered for further transfers)
    worstScen$nests[[i]]$n[worst]=ifelse(worstScen$nests[[i]]$n[worst]<=transfer,0,worstScen$nests[[i]]$n[worst]-transfer)

    #Move TRANSFER foragers out of (and into) worst (and best) cells in BESTNESTS
    bestScen$nests[[i]]$n[best]=bestScen$nests[[i]]$n[best]+transfer #Adds new foragers to best cell
    bestScen$nests[[i]]$n[worst]=bestScen$nests[[i]]$n[worst]-transfer #Subtracts foragers from worst cell

    #Step 3: Update values for worstNests and bestNests

    use=worst|best #Location of worst and best cells

    #For bestNests
    temp=lapply(which(use),optimLoadCurr,scenario=bestScen) #Run optimLoadCurr on all changed cells
    bestScen$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
    for(u in 1:length(temp)){ #For each cell processed
      for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
        bestScen$nests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
        bestScen$nests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns currency
      }
    }

    #For worstNests:
    temp=lapply(which(use),optimLoadCurr,scenario=worstScen) #Run optimLoadCurr on all changed cells
    worstScen$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
    for(u in 1:length(temp)){ #For each cell processed
      for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
        worstScen$nests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
        worstScen$nests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns currency
      }
    }
  }

  #Puts scenarios back into scenario set
  scenarioSet$best=bestScen
  scenarioSet$base=baseScen
  if(!sol) scenarioSet$worst=worstScen

  return(scenarioSet)
}
