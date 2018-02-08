#'@title Move foragers from one cell to another, and re-calculate L, S, and curr
#'
#'@description \code{moveForagers} moves foragers from one cell to another, and
#'  re-calculates L, S, and curr.
#'
#'
#'@param scenarioSet Scenario set list, containing best, base, and worst
#'  scenarios (worst is NA for solitary foragers)
#'@param i Which nest should be used? Currently defaults to 1.
#'@param moves Move list from whichMoves, telling which cells to move from and
#'  to.

#'@return Updated scenario set list
#'
#'@examples
#'#Create scenario
#'nests<-list(nest1=list(xloc=1,yloc=1,n=matrix(c(5,10,15),1),
#'                     whatCurr="eff",sol=T,
#'                     eps=0,L_max=59.5,v=7.8,
#'                     beta=0.102,p_i=1,
#'                     h=matrix(rep(1.5,3),1),
#'                     c_f=0.05,c_i=0.0042,
#'                     H=100,d=matrix(c(50,100,150),1),
#'                     L=matrix(rep(59.5,3),1),
#'                     curr=matrix(c(5,10,15),1),
#'                     steps=5,stepNum=1))
#'
#'world<-list(mu=matrix(rep(8.33e-05,3),1),
#'          flDens=matrix(rep(520,3),1),
#'          e=matrix(rep(14.3,3),1),
#'          l=matrix(rep(1,3),1),
#'          f=matrix(rep(0.86,3),1),
#'          cellSize=10,patchLev=FALSE,
#'          S=matrix(rep(1,3),1),
#'          forageType='random')
#'
#'baseScen=list(nests=nests,world=world) #Baseline scenario
# use=c(1,2,3) #Cells to use
# temp=lapply(use,optimLoadCurr,scenario=baseScen) #Get optimL, currency, and S values from baseline
# baseScen$world$S[use]=sapply(temp,function(x) x$S) #Assigns S-values to baseline scenario
# for(u in 1:length(temp)){ #For each cell processed
#   baseScen$nests[[1]][['L']][use[u]]=temp[[u]][['optimL']][[1]] #Assigns L
#   baseScen$nests[[1]][['curr']][use[u]]=temp[[u]][['optimCurr']][[1]] #Assigns currency
# }
# #Create scenario set
# bestScen=makeBest(baseScen,1)
# worstScen=makeWorst(baseScen,1)
# scenSet=list(best=bestScen,base=baseScen,worst=worstScen)

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
    baseScen$nests[[1]]$L[which(worst)]=sapply(temp,function(x) x$optimL) #Assigns L
    baseScen$nests[[1]]$curr[which(worst)]=sapply(temp,function(x) x$optimCurr) #Assigns currency

    #For bestNests:
    use=worst|best #Location of worst and best cells
    temp=lapply(which(use),optimLoadCurr,scenario=bestScen)
    bestScen$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
    bestScen$nests[[1]]$L[which(use)]=sapply(temp,function(x) x$optimL) #Assigns L
    bestScen$nests[[1]]$curr[which(use)]=sapply(temp,function(x) x$optimCurr) #Assigns currency


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
    # worstScen$nests[[i]]$n[best]=worstScen$nests[[i]]$n[best]+transfer
    worstScen$nests[[i]]$n[best]=baseScen$nests[[i]]$n[best]-transfer #Adds new foragers to best cell
    #Subtracts foragers from worst cell - if < transfer, sets foragers to 0 (won't be considered for further transfers)
    worstScen$nests[[i]]$n[worst]=ifelse(worstScen$nests[[i]]$n[worst]>=transfer,worstScen$nests[[i]]$n[worst]-transfer,0)

    #Move TRANSFER foragers out of (and into) worst (and best) cells in BESTNESTS
    #Adds new foragers to best and worst cells
    bestScen$nests[[i]]$n[best]=baseScen$nests[[i]]$n[best]+transfer
    bestScen$nests[[i]]$n[worst]=baseScen$nests[[i]]$n[worst]+transfer

    #Step 3: Update values for worstNests and bestNests

    use=worst|best #Location of worst and best cells

    #For bestNests
    temp=lapply(which(use),optimLoadCurr,scenario=bestScen) #Run optimLoadCurr on all changed cells
    bestScen$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
    bestScen$nests[[1]]$L[which(use)]=sapply(temp,function(x) x$optimL) #Assigns L
    bestScen$nests[[1]]$curr[which(use)]=sapply(temp,function(x) x$optimCurr) #Assigns currency

    #For worstNests:
    temp=lapply(which(use),optimLoadCurr,scenario=worstScen) #Run optimLoadCurr on all changed cells
    worstScen$world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
    worstScen$nests[[1]]$L[which(use)]=sapply(temp,function(x) x$optimL) #Assigns L
    worstScen$nests[[1]]$curr[which(use)]=sapply(temp,function(x) x$optimCurr) #Assigns currency

  }

  #Puts scenarios back into scenario set
  scenarioSet$best=bestScen
  scenarioSet$base=baseScen
  if(!sol) scenarioSet$worst=worstScen

  return(scenarioSet)
}
