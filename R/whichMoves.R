#' Which moves should be made
#'
#' Helper function to decide if any moves should be made, and if so, where
#' should foragers be transfered FROM and TO.
#'
#' @param scenarioSet Scenario set for current simulation.
#'
#' @return List with three elements, \code{move}, \code{from}, and \code{to}.
#' @export
#'
#' @examples
#' #Create scenario
#' nests<-list(xloc=1,yloc=1,n=matrix(c(5,10,15),1),
#'                    whatCurr="eff",sol=T,
#'                    eps=0,L_max=59.5,v=7.8,
#'                    beta=0.102,p_i=1,
#'                       h=matrix(rep(1.5,3)),
#'                       c_f=0.05,c_i=0.0042,
#'                       H=100,d=matrix(c(50,100,150),1),
#'                       L=matrix(rep(59.5,3),1),
#'                       curr=matrix(c(5,10,15),1),
#'                       steps=5,stepNum=1)
#'world<-list(mu=matrix(rep(8.33e-05,3),1),
#'            flDens=matrix(rep(520,3),1),
#'            e=matrix(rep(14.3,3),1),
#'            l=matrix(rep(1,3),1),
#'            f=matrix(rep(0.86,3),1),
#'            cellSize=10,forageType='omniscient',
#'            S=matrix(rep(1,3),1))
#'baseScen=list(nests=nests,world=world)
#'use=c(1,2,3)
#'temp=lapply(use,optimLoadCurr,scenario=baseScen)
#'baseScen$world$S[use]=sapply(temp,function(x) x$S) #Assigns S-values
#'baseScen$nests$L[use]=sapply(temp,function(x) x$optimL) #Assigns L
#'baseScen$nests$curr[use]=sapply(temp,function(x) x$optimCurr) #Assigns currency

#'#Create scenario set
#'bestScen=makeBest(baseScen)
#'worstScen=makeWorst(baseScen)
#'scenSet=list(best=bestScen,base=baseScen,worst=worstScen)
#'(moves=whichMoves(scenSet)) #What moves should foragers make?

whichMoves=function(scenarioSet=NA){ #Return list with best and worst cells

  #Unpack scenario
  bestNests=scenarioSet$best$nests
  nests=scenarioSet$base$nests
  sol=nests$sol #Solitary foraging
  if(!sol) worstNests=scenarioSet$worst$nests
  transfer=nests$steps[nests$stepNum]
  move=F #Should the move be made?

  if(sol) { #If foragers are solitary

    #Step 1: Find cell to move foragers TO (highest potential currency)
    best=which.max(bestNests$curr) #Location of first best cell
    use=matrix(F,nrow(nests$n),ncol(nests$n)) #Location matrix
    use[best]=T
    best=use #Location of best cell

    #Step 2: Find current cells that have worse currency than best cell
    nests$curr[nests$n<transfer]=NA #Set currency of cells with less than TRANSFER foragers to NA
    potential=which(nests$curr<bestNests$curr[best]) #Current cells with currency less than best cell

    #Step 3:
    if(length(potential)>0){ #If there are worse cells
      worst=potential[which.min(nests$curr[potential])] #Find cell that would benefit the most from moving
      use=matrix(F,nrow(nests$n),ncol(nests$n)) #Location matrix
      use[worst]=T
      worst=use #Location of worst cell
      move=T #Move foragers
    } else { #If no cells are worse, it is impossible to do better
      move=F #Don't move foragers
      worst=NA
      best=NA
    }
    # # OLD METHOD
    # #Step 1: find cell to move foragers FROM
    #
    # #Currency in cells where TRANSFER foragers can be subtracted
    # transferable=ifelse(nests$n>=transfer,nests$curr,NA)
    # worst=which(transferable==min(transferable,na.rm=T),arr.ind=T) #Transferable cell with worst currency
    # if(length(worst)>2) worst=worst[1,] #If there are more than 1 worst cells, choose the first
    # use=matrix(F,nrow(nests$n),ncol(nests$n)) #Location matrix
    # use[worst[1],worst[2]]=T
    # worst=use #Location of worst cell
    #
    # #Step 2: Find cell to move foragers TO
    # worstCurr=bestNests$curr[worst] #Currency in worst cell
    # bestNests$curr[worst]=NA #Sets worst cell to NA (temporarily) so it won't be considered
    # best=which(bestNests$curr==max(bestNests$curr,na.rm=T),arr.ind=T) #Location of first best cell
    #
    # if(length(best)>2) best=best[1,] #If there are more than 1 best cells, choose the first
    # bestNests$curr[worst]=worstCurr #Resets currency in worst cell
    # use=matrix(F,nrow(nests$n),ncol(nests$n)) #Location matrix
    # use[best[1],best[2]]=T
    # best=use #Location of best cell
    #
    # #Step 3:
    # #Is moving to the best cell a better deal than staying in the worst cell?
    # #(i.e. all cells are approximately equal in their currency)?
    # if(bestNests$curr[best]-nests$curr[worst]>nests$eps) {
    #   move=T
    # } else {
    #   move=F
    #   worst=NA
    #   best=NA
    # }

  } else { #If foragers are social

    #Step 1: find cell to move foragers FROM

    #Currency intake in Current and -TRANSFER situation.
    diff1=nests$curr*nests$n-worstNests$curr*worstNests$n

    #Changes all NaNs to zero (in case of a -Inf+Inf situation, which arises when comparing foragers in a worthless cell versus fewer foragers in a worthless cell - summed currency difference b/w cells should still be 0)
    diff1[is.nan(diff1)]=0

    #Currency difference in cells where TRANSFER foragers can be subtracted
    transferable=ifelse(nests$n>=transfer,diff1,NA)
    #Location of cell with the least effect of subtracting TRANSFER foragers - Worst cell (best cell to move FROM)
    worst=which(transferable==min(transferable,na.rm=T),arr.ind=T)

    if(nrow(worst)>1) worst=worst[which.max(nests$d[worst]),] #If there are more than 1 worst cells, choose the furthest
    use=matrix(F,nrow(nests$n),ncol(nests$n)) #Location matrix
    use[worst[1],worst[2]]=T
    worst=use #Location of worst cell

    #Step 2: Find cell to move foragers TO
    worstCurr=bestNests$curr[worst] #Saves currency in worst cell
    bestNests$curr[worst]=NA #Sets worst cell to NA so it won't be considered
    #Currency intake in +TRANSFER situation - Current currency intake
    diff2=bestNests$curr*bestNests$n-nests$curr*nests$n
    #Location of cell with the greatest effect of adding TRANSFER foragers - best cell to move TO
    best=which(diff2==max(diff2,na.rm=T),arr.ind=T)
    if(nrow(best)>1) best=best[which.min(nests$d[best]),] #If there are more than 1 best cells, choose the nearest
    bestNests$curr[worst]=worstCurr #Resets currency
    use=matrix(F,nrow(nests$n),ncol(nests$n)) #Location matrix
    use[best[1],best[2]]=T
    best=use #Location of best cell

    #Step 3:
    #Is the improvement caused by moving TRANSFER foragers >0?

    change=diff2[best]-diff1[worst] #Change in currency intake for the colony
    if(change>nests$eps){
      move=T
    } else {
      move=F
      worst=NA
      best=NA
    }
  }
  return(list(move=move,from=worst,to=best))
}
