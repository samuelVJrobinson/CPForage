forageMod=function(world,nests,iterlim=5000,verbose=F,parallel=F){
  #Internal functions
  if(verbose) print('Starting setup...')


  decimalplaces <- function(x) { #Convenience function for finding number of decimal places
    if ((x %% 1) != 0) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }

  #ERROR HANDLING:
  for(i in 1:length(nests)){ #For each nest, check entry values
    if(any(!(names(nests1[[i]]) %in% c("xloc","yloc","n","whatCurr","sol","constants","eps")))){
      stop('Each CPF nest requires the following arguments:\n "xloc","yloc","n","whatCurr","sol","constants","eps"')
    }
    stopifnot(is.numeric(nests[[i]]$xloc),nests[[i]]$xloc>0,decimalplaces(nests[[i]]$xloc)==0) #X-location
    stopifnot(is.numeric(nests[[i]]$yloc),nests[[i]]$yloc>0,decimalplaces(nests[[i]]$yloc)==0) #Y-location
    stopifnot(is.numeric(nests[[i]]$n),nests[[i]]$n>0,decimalplaces(nests[[i]]$n)==0) #Number of foragers
    stopifnot(is.character(nests[[i]]$whatCurr),nests[[i]]$whatCurr=='eff'|nests[[i]]$whatCurr=='rat') #Currency
    stopifnot(is.logical(nests[[i]]$sol)) #Solitary/social
    stopifnot(is.list(nests[[i]]$constants)) #Forager constants
    stopifnot(is.numeric(nests[[i]]$eps),nests[[i]]$eps>0) #Eps term
    if(any(!(names(nests[[i]]$constants) %in% c("L_max","v","beta","p_i","h","c_f","c_i","H")))){
      stop('Forager constants for each CPF nest require the following arguments:\n "L_max" "v" "beta" "p_i" "h" "c_f" "c_i" "H"')
    }
    if(any(sapply(nests[[i]]$constants,function(x) any(!c(is.numeric(x),x>0))))){
      stop('Forager constants must be numeric, and >0')
    }
  }
  #Check values for world
  if(any(!(names(world) %in% c("mu","flDens","e","l","f","cellSize","patchLev")))){
    stop('Each world requires the following arguments:\n "mu" "flDens" "e" "l" "f" "cellSize" "patchLev"')
  }
  stopifnot(!any(!sapply(world[c('mu','flDens','e','l','f')],is.matrix)), #Are matrices appropriate?
            !any(!sapply(world[c('mu','flDens','e','l','f')],is.numeric)),
            !any(!sapply(world[c('mu','flDens','e','l','f')],function(x) min(x)>=0)))
  stopifnot(length(unique(sapply(world[c('mu','flDens','e','l','f')],ncol)))==1, #Dimension checking
            length(unique(sapply(world[c('mu','flDens','e','l','f')],nrow)))==1)
  stopifnot(is.numeric(world$cellSize),world$cellSize>0,is.logical(world$patchLev))

  #SETUP
  #Add matrix of competition values to the world
  world$S=matrix(1,nrow=nrow(world$mu),ncol=ncol(world$mu))
  world$S[world$mu==0]=NA

  for(i in 1:length(nests)){ #For each nest, do set-up
    nests[[i]]=c(nests[[i]],nests[[i]]$constants)
    nests[[i]]$constants=NULL
    nforagers=nests[[i]]$n #Gets number of foragers for the nest
    #Creates empty vector of forager numbers at each location in the world
    nests[[i]]$n=matrix(0,nrow(world[[1]]),ncol(world[[1]]))
    nests[[i]]$n[nests[[i]]$yloc,nests[[i]]$xloc]=nforagers #Places foragers next to nest
    #Calculates absolute distance of each cell from nest
    nests[[i]]$d=sqrt((((row(world[[1]])-nests[[i]]$yloc)*world$cellSize)^2)+(((col(world[[1]])-nests[[i]]$xloc)*world$cellSize)^2))
    #Empty matrix of loading rate, Load, and currency to be maximized (Rate or Eff)
    nests[[i]]$L=matrix(0,nrow(world[[1]]),ncol(world[[1]])) #Load size (L)
    nests[[i]]$L[nests[[i]]$yloc,nests[[i]]$xloc]=nests[[i]]$L_max #Sets Load size to maximum (initially)
    nests[[i]]$curr=matrix(0,nrow(world[[1]]),ncol(world[[1]])) #Currency
    htemp=nests[[i]]$h
    nests[[i]]$h=matrix(NA,nrow(world[[1]]),ncol(world[[1]]))
    nests[[i]]$h[world$mu>0]=htemp #Handling time - TEMPORARY: FOR FUTURE SCENARIOS, SHOULD BE LOOKED UP IN A TABLE
    #Number of foragers to move at each time step
    if(!is.null(nests[[i]]$steps)){ #If the step number is defined by the user
      #Stops execution if steps are not defined properly
      if(!is.numeric(nests[[i]]$steps)) stop('Step number not properly defined')
      nests[[i]]$steps=sort(nests[[i]]$steps,T)
    } else { #Automatic determination of step size
      #Calculate max number of foragers to move during a single time step (avoiding "reflection" problem, where large number of foragers are assigned to far end of world simply because of depletion effect)
      #NOTE: since implementation of scalFun, there is seems to be no "reflection" problem, so using 5% of nForagers as maxN
      maxN=nforagers*0.05
      nests[[i]]$steps=round(nforagers*1/(10^(seq(1,10,0.5)))) #Initial distribution
      nests[[i]]$steps=c(nests[[i]]$steps[nests[[i]]$steps>1],1) #Gets rid of numbers less than 2, and adds a 1 to the end
      #Cuts off anything step size above maxN, making maxN the largest possible step
      nests[[i]]$steps=c(floor(maxN),nests[[i]]$steps[nests[[i]]$steps<maxN])
      nests[[i]]$steps=sort(unique(nests[[i]]$steps),T) #Descending unique values only
      rm(maxN) #Cleanup
    }
    nests[[i]]$stepNum=1 #Starting point for the steps
  }

  if(parallel){
    if(verbose) cat('Loading parallel processing libraries...')
    require(doSNOW) #Parallel processing
    cluster=makeCluster(4, type = "SOCK") #Set up clusters for parallel processing (4 cores)
    registerDoSNOW(cluster) #Registers clusters
    .Last <- function(){ #"Emergency" function for closing clusters upon unexpected stop
      stopCluster(cluster) #Stops SOCK clusters
      print("forageMod stopped unexpectedly. Closing clusters.")
    }
  } else cluster=NA

  if(verbose) cat('Initializing nests...')

  #TO DO: THIS NEEDS TO BE MODIFIED SO THAT CORRECT L, CURR, S VALUES FROM TEMP ARE ASSIGNED IN THE PROPER NEST (IN MULTI-NEST SCENARIO). CHECK MAKEBEST, MAKEWORST, AND THE REST OF THE MAIN WHILE LOOP TO MAKE SURE THIS IS GOING ON PROPERLY THERE TOO.
  for(i in 1:length(nests)){ #Calculates initial loading rate, load size, and currency for each nest
    #Calculates currency for each cell for each nest
    occupied=nests[[i]]$n>0
    temp=optimLoadCurr(occupied,nests,world) #Optimizes Load and Rate for occupied cells in nest i
    nests[[i]]$L[occupied]=temp$optimL
    nests[[i]]$curr[occupied]=temp$optimCurr
    world$S[occupied]=temp$S
  }

  #worstNests: represents world -transfer foragers in each cell
  #NOTE: CURRENTLY THIS WORKS ONLY FOR SINGLE-NEST SITUATION.
  if(sum(sapply(nests, function(x) !x$sol))>0){ #If at least one nest uses social foraging
    temp=makeWorst(nests,world,whichNest=1,cluster=cluster)
    worstNests=temp$worstNests
    worstWorld=temp$worstWorld
  }

  #NOTE: CURRENTLY THIS WORKS ONLY FOR SINGLE-NEST SITUATION. SHOULD BE LOOPED OVER ALL COMBINATIONS OF 'MOVES' IN MULTI-NEST SITUATION. IT WOULD PROBABLY HELP IF WORLD + NESTS WAS GROUPED TOGETHER IN LISTS ('SCENARIOS')
  #bestNests: represents world +transfer foragers in each cell
  temp=makeBest(nests,world,whichNest=1,parallel=parallel,cluster=cluster)
  bestNests=temp$bestNests
  bestWorld=temp$bestWorld

  #Index of which nests are "done" (can't improve distribution any further)
  done=rep(F,length(nests))
  nitt=1 #Number of iterations (debugging to see when it should be "cut off")

  if(verbose) print(paste('Simulation started at',Sys.time()))
  #Main loop
  while(sum(!done)>0){
    if(verbose && (nitt %% 10)==0) print(paste('Iteration',nitt)) #Prints every 10 iterations
    done=rep(F,length(nests)) #Resets each time that at least one is "not done"
    for(i in 1:length(nests)){  #For each nest
      transfer=nests[[i]]$steps[nests[[i]]$stepNum] #Number of foragers to transfer
      if(nests[[i]]$sol){ #If optimization is being done for solitary foragers

        ###Figure out where worst and best cells are ###

        worst=which(ifelse(nests[[i]]$n>transfer-1,nests[[i]]$curr,NA)==min(ifelse(nests[[i]]$n>transfer-1,nests[[i]]$curr,NA),na.rm=T),arr.ind=T) #Worst cell
        if(length(worst)>2) worst=worst[1,] #If there are more than 1 worst cells, choose the first
        use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
        use[worst[1],worst[2]]=T
        worst=use #Location of worst cell

        worstCurr=bestNests[[i]]$curr[worst] #Currency in worst cell
        bestNests[[i]]$curr[worst]=NA #Sets worst cell to NA (temporarily) so it won't be considered
        best=which(bestNests[[i]]$curr==max(bestNests[[i]]$curr,na.rm=T),arr.ind=T) #Location of first best cell
        if(length(best)>2) best=best[1,] #If there are more than 1 best cells, choose the first
        bestNests[[i]]$curr[worst]=worstCurr #Resets currency in worst cell
        use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
        use[best[1],best[2]]=T
        best=use #Location of best cell

        #Is moving to the best cell a better deal than staying in the worst cell?
        #(i.e. all cells are approximately equal in their currency)?
        if(bestNests[[i]]$curr[best]-nests[[i]]$curr[worst]>nests[[i]]$eps){

          #Move TRANSFER foragers to the best cell
          nests[[i]]$n[best]=bestNests[[i]]$n[best]
          #Overwrite optimal Load and Currency from the best cell
          nests[[i]]$L[best]=bestNests[[i]]$L[best] #Load size
          nests[[i]]$curr[best]=bestNests[[i]]$curr[best] #Currency
          world$S[best]=bestWorld$S[best] #S-value

          #Move TRANSFER foragers out of worst cell in NESTS
          nests[[i]]$n[worst]=nests[[i]]$n[worst]-transfer
          if(nests[[i]]$n[worst]==0 & sum(sapply(nests,function(x) x$n[worst]))==0){ #If worst cell is now empty...
            nests[[i]]$L[worst]=nests[[i]]$curr[worst]=0 #Set Load and currency to 0
            world$S[worst]=1 #Set S to 1 (no competitors)
          } else {
            #Calculate new optimal load and currency for worst cell in NESTS
            temp=optimLoadCurr(u=worst,nests,world)
            world$S[worst]=temp$S #S-value
            for(name in names(temp$optimCurr)){
              nests[[name]][['L']][worst]=temp[['optimL']][[name]] #Assigns L
              nests[[name]][['curr']][worst]=temp[['optimCurr']][[name]] #Assigns currency
            }
          }
          #Move TRANSFER foragers out of (and into) worst (and best) cells in BESTNESTS
          bestNests[[i]]$n[best]=bestNests[[i]]$n[best]+transfer #Adds new foragers to best cell
          bestNests[[i]]$n[worst]=bestNests[[i]]$n[worst]-transfer #Subtracts foragers from worst cell
          #Calculate new optimal load and currency for worst and best cells in BESTNESTS
          use=worst|best #Location of worst and best cells
          temp=lapply(which(use),optimLoadCurr,nests=bestNests,world=bestWorld)
          world$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
          for(u in 1:length(temp)){ #For each cell processed
            for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
              bestNests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
              bestNests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns currency
            }
          }
          #End of IF statement
        } else if(transfer>1) { #If transfer number is >1 (i.e. not at the end of the list)
          if(verbose) print(paste('Finished pass',nests[[i]]$stepNum,'of',length(nests[[i]]$steps),'for nest',i,'. Starting pass ',nests[[i]]$stepNum+1,'...'))
          nests[[i]]$stepNum=nests[[i]]$stepNum+1 #Increments step number
          #Creates new "bestNest" and "bestWorld" (for nest i)
          temp=makeBest(nests,world,whichNest=i,parallel=parallel,cluster=cluster)
          bestNests=temp$bestNests
          bestWorld=temp$bestWorld
        } else {done[i]=T} #If transfer is 1, and there's no better deal, distribution has converged
        #SOCIAL FORAGERS
      } else { #If optimization is being done for social foragers, colony rate/efficiency is maximized

        #Step 1: find cell to move foragers FROM

        #THIS IS FAILING BECAUSE CURRENCY CALCULATIONS IN EMPTY CELLS RETURN -INF FOR RATE. -INF+INF = Nan
        #Currency intake in Current and -TRANSFER situation
        diff1=nests[[i]]$curr*nests[[i]]$n-worstNests[[i]]$curr*worstNests[[i]]$n
        #Changes all NaNs to zero (in case of a -Inf+Inf situation, which arises when comparing foragers in a worthless cell versus fewer foragers in a worthless cell - summed currency difference b/w cells should still be 0)
        diff1[is.nan(diff1)]=0
        #Location of cell with the least effect of subtracting TRANSFER foragers - Worst cell
        worst=which(ifelse(nests[[i]]$n>transfer-1,diff1,NA)==min(ifelse(nests[[i]]$n>transfer-1,diff1,NA),na.rm=T),
                    arr.ind=T)
        if(length(worst)>2) worst=worst[1,] #If there are more than 1 worst cells, choose the first
        use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
        use[worst[1],worst[2]]=T
        worst=use #Location of worst cell

        #Step 2: Find cell to move foragers TO
        worstCurr=bestNests[[i]]$curr[worst] #Saves currency in worst cell
        bestNests[[i]]$curr[worst]=NA #Sets worst cell to NA so it won't be considered
        #Currency intake in +TRANSFER situation - Current currency intake
        diff2=bestNests[[i]]$curr*bestNests[[i]]$n-nests[[i]]$curr*nests[[i]]$n
        best=which(diff2==max(diff2,na.rm=T),arr.ind=T) #Location of cell with the greatest effect of adding TRANSFER foragers
        if(length(best)>2) best=best[1,] #If there are more than 1 best cells, choose the first
        bestNests[[i]]$curr[worst]=worstCurr #Resets currency
        use=matrix(F,nrow(nests[[i]]$n),ncol(nests[[i]]$n)) #Location matrix
        use[best[1],best[2]]=T
        best=use #Location of best cell

        #THIS IS RETURNING MISSING ARGUMENTS, MEANING THAT CIRCUMSTANCES WHERE NO >TRANSFER CELLS WERE PASSED TO IT.

        #Step 3: Is the improvement caused by moving TRANSFER foragers >0?
        change=diff2[best]-diff1[worst] #Change in currency intake for the colony
        if(change>nests[[i]]$eps){
          #Move TRANSFER foragers from worst cell
          nests[[i]]$n[worst]=worstNests[[i]]$n[worst] #Overwrite forager number for worst cell
          nests[[i]]$L[worst]=worstNests[[i]]$L[worst] #Load size
          nests[[i]]$curr[worst]=worstNests[[i]]$curr[worst] #Currency

          #Move TRANSFER foragers to best cell
          nests[[i]]$n[best]=bestNests[[i]]$n[best] #Overwrite forager number for best cell
          nests[[i]]$L[best]=bestNests[[i]]$L[best] #Load size
          nests[[i]]$curr[best]=bestNests[[i]]$curr[best] #Currency

          #Move TRANSFER foragers out of (and into) worst (and best) cells in BESTNESTS
          bestNests[[i]]$n[best]=bestNests[[i]]$n[best]+transfer #Adds new foragers to best cell
          bestNests[[i]]$n[worst]=bestNests[[i]]$n[worst]-transfer #Subtracts foragers from worst cell
          #Calculate new optimal load and currency for worst and best cells in BESTNESTS
          use=worst|best #Location of worst and best cells
          temp=lapply(which(use),optimLoadCurr,nests=bestNests,world=bestWorld)
          bestWorld$S[which(use)]=sapply(temp,function(x) x$S) #Assigns S-values
          for(u in 1:length(temp)){ #For each cell processed
            for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
              bestNests[[name]][['L']][which(use)[u]]=temp[[u]][['optimL']][[name]] #Assigns L
              bestNests[[name]][['curr']][which(use)[u]]=temp[[u]][['optimCurr']][[name]] #Assigns currency
            }
          }

          #Move TRANSFER foragers out of (and into) worst (and best) cells in WORSTNESTS
          worstNests[[i]]$n[best]=worstNests[[i]]$n[best]+transfer #Adds new foragers to best cell
          worstNests[[i]]$n[worst]=worstNests[[i]]$n[worst]-transfer #Subtracts foragers from worst cell
          #If worst cell is now empty...
          if(worstNests[[i]]$n[worst]==0 & sum(sapply(worstNests,function(x) x$n[worst]))==0){
            worstNests[[i]]$L[worst]=worstNests[[i]]$curr[worst]=0 #Set Load and currency to 0
            worstWorld$S[worst]=1 #Set S to 1 (no competitors)
          } else {
            #Calculate new optimal load and currency for worst cell in worstNests
            temp=optimLoadCurr(u=worst,worstNests,worstWorld)
            worstWorld$S[worst]=temp$S #S-value
            for(name in names(temp$optimCurr)){
              worstNests[[name]][['L']][worst]=temp[['optimL']][[name]] #Assigns L
              worstNests[[name]][['curr']][worst]=temp[['optimCurr']][[name]] #Assigns currency
            }
          }
          #By definition, best cell cannot be completely empty, therefore:
          #Calculate new optimal load and currency for best cell in worstNests
          temp=optimLoadCurr(u=best,worstNests,worstWorld)
          worstWorld$S[best]=temp$S #S-value
          for(name in names(temp$optimCurr)){
            worstNests[[name]][['L']][best]=temp[['optimL']][[name]] #Assigns L
            worstNests[[name]][['curr']][best]=temp[['optimCurr']][[name]] #Assigns currency
          }

        } else if(transfer>1) { #If transfer number is >1 (i.e. not at the end of the list)
          if(verbose) print(paste('Finished pass',nests[[i]]$stepNum,'of',length(nests[[i]]$steps),'for nest',i))

          nests[[i]]$stepNum=nests[[i]]$stepNum+1 #Increments step number

          #Creates new "bestNest" and "bestWorld" (for nest i)
          temp=makeBest(nests,world,whichNest=i,parallel=parallel,cluster=cluster)
          bestNests=temp$bestNests
          bestWorld=temp$bestWorld

          #Creates new "bestNest" and "bestWorld" (for nest i)
          temp=makeWorst(nests,world,whichNest=i,parallel=parallel,cluster=cluster)
          worstNests=temp$worstNests
          worstWorld=temp$worstWorld

        } else {done[i]=T} #If transfer is 1, and there's no better deal, distribution has converged for this next
      } #End of IF(sol/soc)
    } #End of FOR loop
    #If all nests are "done", loop should exit
    nitt=nitt+1 #Increment counter
    if(nitt==iterlim) {
      if(verbose) print('Iteration limit reached')
      break
    }
  } #End of WHILE loop

  #Calculate patch residence times per forager (in seconds)
  for(i in 1:length(nests)){
    nests[[i]]$loadingTime=with(nests[[i]],ifelse(n>0,L*(h+world$S*world$l*p_i+world$f)/world$S*world$l,NA))
    nests[[i]]$travelTime=with(nests[[i]],ifelse(n>0,(d*(2-beta*(L/L_max)))/(v*(1-beta*(L/L_max))),NA))
    nests[[i]]$boutLength=nests[[i]]$loadingTime+nests[[i]]$travelTime+nests[[i]]$H #Time for 1 complete foraging bout
  }
  if(parallel) stopCluster(cluster) #Stops SOCK clusters
  if(verbose) print(paste('Simulation ended at',Sys.time(),'Final number of iterations = ',nitt,'.'))
  return(list(world=world,nests=nests)) #Returns world and nests in a list
}
