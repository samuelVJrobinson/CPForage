whichMoves=function(scenarioSet=NA,i=NA){ #Using nest i in scenarioSet, return list with best and worst cells
  if(length(scenarioSet)==1 & is.na(scenarioSet[1])){
    stop('Nests not specified')
  } else if(is.na(i)){
    stop('Which nest not specified')
  }

  #Unpack scenario
  bestNests=scenarioSet$best$nests
  nests=scenarioSet$base$nests
  sol=nests[[i]]$sol #Solitary foraging
  if(!sol) worstNests=scenarioSet$worst$nests
  transfer=nests[[i]]$steps[nests[[i]]$stepNum]
  move=F #Should the move be made?

  if(sol) { #If foragers are solitary

    #Step 1: find cell to move foragers FROM

    #Currency in cells where TRANSFER foragers can be subtracted
    transferable=ifelse(nests[[i]]$n>=transfer,nests[[i]]$curr,NA)
    worst=which(transferable==min(transferable,na.rm=T),arr.ind=T) #Transferable cell with worst currency
    if(length(worst)>2) worst=worst[1,] #If there are more than 1 worst cells, choose the first
    use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
    use[worst[1],worst[2]]=T
    worst=use #Location of worst cell

    #Step 2: Find cell to move foragers TO
    worstCurr=bestNests[[i]]$curr[worst] #Currency in worst cell
    bestNests[[i]]$curr[worst]=NA #Sets worst cell to NA (temporarily) so it won't be considered
    best=which(bestNests[[i]]$curr==max(bestNests[[i]]$curr,na.rm=T),arr.ind=T) #Location of first best cell

    if(length(best)>2) best=best[1,] #If there are more than 1 best cells, choose the first
    bestNests[[i]]$curr[worst]=worstCurr #Resets currency in worst cell
    use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
    use[best[1],best[2]]=T
    best=use #Location of best cell

    #Step 3:
    #Is moving to the best cell a better deal than staying in the worst cell?
    #(i.e. all cells are approximately equal in their currency)?
    if(bestNests[[i]]$curr[best]-nests[[i]]$curr[worst]>nests[[i]]$eps) {
      move=T
    } else {
      move=F
      worst=NA
      best=NA
    }
  } else { #If foragers are social

    #Step 1: find cell to move foragers FROM

    #Currency intake in Current and -TRANSFER situation.
    diff1=nests[[i]]$curr*nests[[i]]$n-worstNests[[i]]$curr*worstNests[[i]]$n

    #Changes all NaNs to zero (in case of a -Inf+Inf situation, which arises when comparing foragers in a worthless cell versus fewer foragers in a worthless cell - summed currency difference b/w cells should still be 0)
    diff1[is.nan(diff1)]=0

    #Currency difference in cells where TRANSFER foragers can be subtracted
    transferable=ifelse(nests[[i]]$n>=transfer,diff1,NA)
    #Location of cell with the least effect of subtracting TRANSFER foragers - Worst cell (best cell to move FROM)
    worst=which(transferable==min(transferable,na.rm=T),arr.ind=T)

    if(length(worst)>2) worst=worst[1,] #If there are more than 1 worst cells, choose the first
    use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
    use[worst[1],worst[2]]=T
    worst=use #Location of worst cell

    #Step 2: Find cell to move foragers TO
    worstCurr=bestNests[[i]]$curr[worst] #Saves currency in worst cell
    bestNests[[i]]$curr[worst]=NA #Sets worst cell to NA so it won't be considered
    #Currency intake in +TRANSFER situation - Current currency intake
    diff2=bestNests[[i]]$curr*bestNests[[i]]$n-nests[[i]]$curr*nests[[i]]$n
    #Location of cell with the greatest effect of adding TRANSFER foragers - best cell to move TO
    best=which(diff2==max(diff2,na.rm=T),arr.ind=T)
    if(length(best)>2) best=best[1,] #If there are more than 1 best cells, choose the first
    bestNests[[i]]$curr[worst]=worstCurr #Resets currency
    use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
    use[best[1],best[2]]=T
    best=use #Location of best cell

    #Step 3:
    #Is the improvement caused by moving TRANSFER foragers >0?

    #This equation is wrong. Should be
    change=diff2[best]-diff1[worst] #Change in currency intake for the colony
    if(change>nests[[i]]$eps){
      move=T
    } else {
      move=F
      worst=NA
      best=NA
    }
  }
  return(list(move=move,from=worst,to=best))
}
