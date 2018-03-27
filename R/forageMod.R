#'@title Run CPF model
#'
#'@description \code{forageMod} runs central-place foraging model given forager and world
#'data.
#'
#' Function to run central-place foraging (CPF) model based on the ideal-free
#' distribution (IFD). Takes a list of nest parameters, and a list of world
#' parameters, runs the model until convergence, and then returns a list
#' containing a matrix of competitive effects, and a list of matrices of foraging
#' parameters (e.g. time in patch, foraging currency experienced at each cell).
#' Use \code{nests2df()} to convert this to a more readable dataframe. Currently
#' this only works for individual aggregations of CP foragers.
#'
#'@param world World structure. List.
#'@param nests Nests structure. List of lists.
#'@param iterlim Limit to number of iterations. Default = 5000.
#'@param verbose Should function display progress?
#'@param parallel Should parallel processing be used for large tasks?
#'@param ncore Number of SNOW cores to use (if parallel = TRUE).
#'@param parMethod Message passing for parallel processing (see details below).
#'@param tol Tolerance range for optimization function. Default =
#'  .Machine$double.eps^0.25.

#'@return List containing world structure (competition term) and nest structure
#'  (forager distribution)
#'
#'@details \code{parMethod} must be either \code{'SOCK'} (Default) or \code{'MPI'}. Requires \code{doSNOW} and \code{RMPI} (if using MPI) packages.
#'
#'\code{world} should be a named list containing:
#' \itemize{
#' \item \code{mu}: nectar production values (per s); matrix
#' \item \code{e}: energy value of nectar (J/\eqn{\mu}L); matrix
#' \item \code{l}: maximum nectar standing crop (\eqn{\mu}L); matrix
#' \item \code{f}: travel time between flowers (s); matrix
#' \item \code{flDens}: flower count per cell; matrix
#' \item \code{cellSize}: size of a cell (m)
#' \item \code{forageType}: foraging type. See \code{\link{curr}}.
#' }
#'
#' \code{nests} should be a list of lists, each containing:
#' \itemize{
#' \item \code{xloc}: x-location of nest (column number) in world; integer.
#' \item \code{yloc}: y-location of nest (row number) in world; integer.
#' \item \code{n}: number of foragers; integer.
#' \item \code{whatCurr}: '\code{eff}' (efficiency) or '\code{rat}' (rate).
#' \item \code{sol}: solitary foraging; logical.
#' \item \code{constants}: named list of foraging parameters:
#' \itemize{
#' \item \code{L_max}: maximum load (\eqn{\mu}L); numeric.
#' \item \code{v}: maximum flight speed (m/s); numeric.
#' \item \code{beta}: reduction of flight speed with load (m/s*\eqn{\mu}L); numeric.
#' \item \code{p_i}: rate of nectar uptake ("licking speed", \eqn{\mu}L/s); numeric.
#' \item \code{h}: handling time before draining a flower (s); numeric.
#' \item \code{c_f}: energetic cost of flight (J/s); numeric.
#' \item \code{c_i}: energetic cost of non-flight (J/s); numeric.
#' \item \code{H}: time spent in hive/aggregation (s); numeric.
#' }
#' \item \code{eps}: accuracy to use for optimization; numeric.
#' }
#'
#'
#'@examples
#'
#'#Create test world for run
#'nu_i<-0.3/3600 #Nectar production/hr for a single flower
#'flDens<-520 #Flower density/m2
#'e_i<-14.35 #Energetic value/unit
#'l_i<-1 #Canola standing crop (1uL)
#'f_i<-0.86 #Inter-flower flight time
#
#'#World structure
#'cellSize<-10 #10m cells (100m^2)
#'worldSize<-120 #120x120m field (100x100m field with 10m buffer zone worth nothing)
#'world1<-list(mu=matrix(0,120,120),
#'            flDens=matrix(0,120,120),
#'            e=matrix(0,120,120),
#'            l=matrix(0,120,120),
#'            f=matrix(0,120,120),
#'            cellSize=cellSize) #Empty world
#'world1$mu[c(2:11),c(2:11)]<-nu_i #Per-flower nectar production in
#'canola-filled cells
#'world1$flDens[c(2:11),c(2:11)]<-flDens*cellSize^2 #Flower number per cell
#'world1$e[c(2:11),c(2:11)]<-e_i #Energy production in canola-filled cells
#'world1$l[c(2:11),c(2:11)]<-l_i #Standing crop in cells with no competition
#'world1$f[c(2:11),c(2:11)]<-f_i #Inter-flower flight time world1$patchLev=F
#'world1$forageType <- 'omniscient' #Foraging style for flowers within patch
#'
#'#Constants for foragers
#'honeybeeConstants<-list(L_max=59.5, #Max load capacity (uL) - Schmid-Hempel (1987)
#'                       v=7.8, #Velocity (m/s) - Unloaded flight speed (Wenner 1963)
#'                       beta=0.102, #Proportion reduction in completely loaded
#'                       flight speed (1-v/v_l)
#'                       p_i=1, # Max loading rate (uL/s)
#'                       h=1.5, #Handling time per flower (s)
#'                       c_f=0.05, #Unloaded flight energetic cost (J/s) (Dukas
#'                       and Edelstein Keshet 1998)
#'                       c_i=0.0042, #Cost of non-flying activity
#'                       H=100) #Time spent in the hive (s) (Seeley 1986 found
#'                       100s and 70s for high and low intake rates)
#'
#'#Nest structure (social rate maximizers)
#'nests1<-list(nest1=list(xloc=1,yloc=1,n=1000,whatCurr='rat',sol=F,constants=honeybeeConstants))
#'
#'#Run model
#'testOutput1<-forageMod(world1,nests1,2000,verbose=F,parallel=T)

forageMod=function(world,nests,iterlim=5000,verbose=F,parallel=F,ncore=4,parMethod='SOCK',tol=.Machine$double.eps^0.25){
  #Internal functions
  if(verbose) print('Starting setup...')
  decimalplaces <- function(x) { #Convenience function for finding number of decimal places
    if ((x %% 1) != 0) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  #INPUT CHECKING:
  if(length(nests)>1){
    stop("forageMod can't do multi-nest simulations (yet).")
  }
  if(is.null(names(nests))){ #Checks for NULL or empty nest names
    if(length(nests)==1) names(nests)='nest1' else names(nests)=paste(rep('nest',length(nests)),1:length(nests),sep='')
    warning('Nests unnamed. Providing names.')
  } else if(any(nchar(names(nests))==0)) {
    #Provides a new name
    names(nests)=ifelse(nchar(names(nests))==0,paste('nests',which(nchar(names(nests))==0),sep=''),names(nests))
    warning('Some nests are unnamed. Providing name for unnamed nests.')
  }
  if(length(names(nests))>1 & any(duplicated(names(nests)))) {
    names(nests)[duplicated(names(nests))]=
      replicate(sum(duplicated(names(nests))),paste(cbind(sample(letters,8)),collapse=''))
    warning('Duplicate nest names. Providing random names to duplicate nests.')
  }

  for(i in 1:length(nests)){ #For each nest, check entry values
    if(any(!(c("xloc","yloc","n","whatCurr","sol","constants","eps") %in% names(nests[[i]])))){
      stop('Each CPF nest requires the following arguments:\n "xloc","yloc","n","whatCurr","sol","constants","eps"')
    stopifnot(is.numeric(nests[[i]]$xloc),nests[[i]]$xloc>0,decimalplaces(nests[[i]]$xloc)==0) #X-location
    }
    stopifnot(is.numeric(nests[[i]]$yloc),nests[[i]]$yloc>0,decimalplaces(nests[[i]]$yloc)==0) #Y-location
    stopifnot(is.numeric(nests[[i]]$n),nests[[i]]$n>0,decimalplaces(nests[[i]]$n)==0) #Number of foragers
    stopifnot(is.character(nests[[i]]$whatCurr),nests[[i]]$whatCurr=='eff'|nests[[i]]$whatCurr=='rat') #Currency
    stopifnot(is.logical(nests[[i]]$sol)) #Solitary/social
    stopifnot(is.list(nests[[i]]$constants)) #Forager constants
    stopifnot(is.numeric(nests[[i]]$eps),nests[[i]]$eps>=0) #Eps term
    if(any(!(names(nests[[i]]$constants) %in% c("L_max","v","beta","p_i","h","c_f","c_i","H")))){
      stop('Forager constants for each CPF nest require the following arguments:\n "L_max" "v" "beta" "p_i" "h" "c_f" "c_i" "H"')
    }
    if(any(sapply(nests[[i]]$constants,function(x) any(!c(is.numeric(x),x>0))))){
      stop('Forager constants must be numeric, and >0')
    }
  }
  #Check values for world
  if(any(!(c("mu","flDens","e","l","f","cellSize","forageType") %in% names(world)))){
    stop('Each world requires the following arguments:\n "mu" "flDens" "e" "l" "f" "cellSize" "forageType"')
  }
  #Are matrices appropriately defined?
  stopifnot(!any(!sapply(world[c('mu','flDens','e','l','f')],is.matrix)),
            !any(!sapply(world[c('mu','flDens','e','l','f')],is.numeric)),
            !any(!sapply(world[c('mu','flDens','e','l','f')],function(x) min(x)>=0)))
  stopifnot(length(unique(sapply(world[c('mu','flDens','e','l','f')],ncol)))==1, #Dimension checking
            length(unique(sapply(world[c('mu','flDens','e','l','f')],nrow)))==1)
  stopifnot(is.numeric(world$cellSize),world$cellSize>0,is.character(world$forageType))

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
    nests[[i]]$d=nests[[i]]$d+min(nests[[i]]$d[nests[[i]]$d!=0])/2 #Adds half minimum distance (prevents 0 flying distance)
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
      fakeNests=nests #Copy of nests
      fakeNests[[i]]$n=0 #No foragers
      fakeNests[[i]]$h=min(fakeNests[[i]]$h,na.rm=T) #minimum handling time
      fakeNests[[i]]$d=min(fakeNests[[i]]$d[fakeNests[[i]]$d!=0],na.rm=T) #minimum nonzero distance
      fakeNests[[i]]$L=fakeNests[[i]]$L_max #L_max
      fakeNests[[i]]$curr=0
      fakeWorld=world #Copy of world
      richest=which.max(fakeWorld$mu*fakeWorld$flDens*fakeWorld$e)
      fakeWorld$mu=fakeWorld$mu[richest] #Mu from richest patch
      fakeWorld$flDens=fakeWorld$flDens[richest] #flDens from richest patch
      fakeWorld$e=fakeWorld$e[richest] #e ...
      fakeWorld$l=fakeWorld$l[richest] #l
      fakeWorld$f=fakeWorld$f[richest] #f
      fakeWorld$S=1 #No competition (initially)
      #Function to return currency in a cell given n foragers
      maxNfun=function(n,fakeNests,fakeWorld,i,eps){
        fakeNests[[i]]$n=n
        temp=optimLoadCurr(1,list(nests=fakeNests,world=fakeWorld))
        return(abs(temp$S-eps)) #Return S-value at given n
      }
      #Number of foragers where S goes to Smin (essentially upper limit to step size)
      #Simulation indicates that Smin should optimally be around 0.18 for multi-core, 0.3 for serial runs.
      maxn=optimize(maxNfun,interval=c(0,max(nests[[i]]$n)),fakeNests=fakeNests,
                    fakeWorld=fakeWorld,i=i,eps=0.5)$min
      #Simulation indicates that phi should be around 4.3, Smin around 0.7
      nests[[i]]$steps=round(nforagers*1/(10^(seq(1,10,4.3)))) #Initial distribution
      nests[[i]]$steps=c(nests[[i]]$steps[nests[[i]]$steps>1],1) #Gets rid of numbers less than 2, and adds a 1 to the end
      #Cuts off anything step size above maxN, making maxN the largest possible step
      nests[[i]]$steps=c(floor(maxn),nests[[i]]$steps[nests[[i]]$steps<maxn])
      nests[[i]]$steps=sort(unique(nests[[i]]$steps),T) #Descending unique values only
      rm(maxn,maxNfun) #Cleanup
    }
    nests[[i]]$stepNum=1 #Starting point for the steps
  }

  if(parallel){
    if(verbose) cat('Loading parallel processing libraries...')
    require(doSNOW) #Parallel processing

    #Set up ncore number of clusters for parallel processing
    if(parMethod=='MPI'){
      require(Rmpi)
      ncore <- mpi.universe.size() #Number of cores available to spawn processes
      sprintf("TEST mpi.universe.size() =  %i", ncore) #Number of MPI processes
      cluster <- makeCluster(ncore-1, type = "MPI") #Set up clusters for parallel processing (# cores - 1 = # slave cores)
    } else {
      cluster <- makeCluster(ncore, type = "SOCK") #Create SOCK clusters
    }

    registerDoSNOW(cluster) #Registers clusters

    #"Emergency" function for closing clusters upon unexpected stop
    .Last <- function(){
      stopCluster(cluster) #Stops SOCK clusters
      if(parMethod=='MPI') mpi.exit()
      print("forageMod stopped unexpectedly. Closing clusters.")
    }
  } else cluster=NA

  if(verbose) cat('Initializing nests...')

  for(i in 1:length(nests)){ #Calculates initial loading rate, load size, and currency for each nest
    #Calculates currency for each cell for each nest
    occupied=nests[[i]]$n>0
    temp=optimLoadCurr(which(occupied),list(nests=nests,world=world)) #Optimizes Load and Rate for occupied cells in nest i
    world$S[occupied]=temp$S #S-value
    nests[[i]][['L']][occupied]=temp[['optimL']][[i]] #Assigns L
    nests[[i]][['curr']][occupied]=temp[['optimCurr']][[i]] #Assigns currency
  }

  #Idea:
  #Create BASE scenario: list containing a world and set of nests-represents "current situation"
  base=list(nests=nests,world=world)

  #Depending on competitive framework:

  #1: SIMPLE CASE - 1 nest
  #create 1 (or 2) extra scenarios for solitary (or social) foragers.
  #Add TRANSFER foragers to 1st scenario, rerun optimLoadCurr to get S/L/currency values (similar to "bestNests")
  #If nest is social:
    #Subtract TRANSFER foragers to 2nd scenario, rerun optimLoadCurr to get S/L/currency values (similar to "worstNests")
  #Compare currency values b/w 1st and BASE scenario (or 1st, BASE, and 2nd scenario)
  #If currency change is better than EPS
  #Change forager number in BEST and WORST cells in all 3 scenarios
  #Run optimLoadCurr for candidate cells in all 3 scenarios
  #Loop

  #2: COMPLEX CASE - 2 or more nests - ignoring other nest behaviour
  #Requires (n*2 + 1) sets of calculations for each iteration

  #Create nests:
  #for(i = nests)
  #Create 1 (or 2) extra scenarios
  #Scenario1 = Add TRANSFER foragers to i in BASE,  rerun optimLoadCurr
  #If nest i is social:
    #Scenario2 = Subtract TRANSFER foragers from i in BASE, rerun optimLoadCurr
  #end

  #Compare nests:
  #for(i = nests)
  #Compare singular currency values b/w 1st, (2nd) and BASE scenario
  #end


  #3: COMPLEX CASE - 2 or more nests - marginalizing across other nest behaviour
  #Requires (n*3^(n-1)) sets of calculations for each iteration
  #Create nests:
  #for(i = nests)
  #Create 2 (or 3) extra scenario SETS
  #ScenarioSet0 = Use i in BASE
  #= Create combination of ADD/ZERO(/SUBTRACT) TRANSFER for all other nests
  #= Note: 1 of these scenarios can just be copied from BASE
  #= Run all scenarios through optimLoadCurr
  #ScenarioSet1 = Add TRANSFER foragers to i in BASE
  #= Create combination of ADD/ZERO(/SUBTRACT) TRANSFER for all other nests
  #= Run all scenarios through optimLoadCurr
  #If nest i is social:
    #Scenario2 = Subtract TRANSFER foragers from i in BASE
    #= Create combination of ADD/ZERO(/SUBTRACT) TRANSFER for all other nests (requires 3^(n-1) scenarios for all combinations)
    #= Run all scenarios through optimLoadCurr
  #end

  #Compare nests:
  #for(i = nests)
    #Compare multiple currency values b/w 'best', ('worst') and 'zero' scenarios.
    #Fetch single value for each cell in the entire "stack" of scenarios. Options:
      #Summed currencies
      #Mean value of currencies
      #Median value of currrencies
      #Maximum currency (best-case scenario)
      #Minimum currency (worst-case scenario)
      #Other extensions: different functions for each scenario set? (e.g. maximum from 'zero', but minimum from 'best')
    #For solitary: choose max(best), and min(zero) for best and worst cells
    #For social: choose max(diff1), and min(diff2) for best and worst cells
  #Etc...
  #end

  #REQUIREMENTS:
  #1) Set up model to deal with scenario=list(world,nests)
    #a) Modify optimLoadCurr to deal with & return scenarios
    #b) Any other functions underneath optimLoadCurr?? (Don't think they need modification)
  #2) Set up model preamble to deal with 'ignore' or 'marginalize' behaviour
    #a) BASE scenario should be in environment
    #b) If('ignore') each nest should have a list with 1 (or 2) scenarios, "best" (and "worst")
    #c) If('marginalize') each nest should have a list with 1 (or 2) scenario SETS, "best", "zero" (and "worst").
    #"best" should contain all scenarios


  #NOTE: CURRENTLY THIS WORKS ONLY FOR SINGLE-NEST SITUATION. SHOULD BE LOOPED OVER ALL COMBINATIONS OF 'MOVES' IN MULTI-NEST SITUATION. ALTERNATIVELY, IT COULD JUST USE MARGINAL SITUATION (WHAT MOVE | OTHERS), BUT THIS MAY LEAD TO STABLE OSCILLATIONS, OR TURN-BASED ADVANTAGES

  #worstNests: represents world -transfer foragers in each cell
  if(sum(sapply(base$nests, function(x) !x$sol))>0){ #If at least one nest uses social foraging
    worst=makeWorst(base,whichNest=1,cluster=cluster)
  } else {
    worst=NA
  }

  #bestNests: represents world +transfer foragers in each cell
  best=makeBest(base,whichNest=1,parallel=parallel,cluster=cluster)

  #Groups scenarios (best, base, worst) into scenario set for nest1
  nestSet=list(best=best,base=base,worst=worst)

  nNests=length(nestSet$base$nests) #Number of nests involved

  #Index of which nests are "done" (can't improve distribution any further)
  done=rep(F,nNests)
  nitt=1 #Number of iterations (debugging to see when it should be "cut off")
  startTime <- Sys.time() #Starting time
  if(verbose) print(paste('Simulation started at',startTime))
  #Main loop
  while(sum(!done)>0){
    if(verbose && (nitt %% 10)==0) print(paste('Iteration',nitt)) #Prints every 10 iterations
    done=rep(F,length(nestSet$base$nests)) #Resets each time that at least one is "not done"
    for(i in 1:nNests){  #For each nest
      transfer=with(nestSet$base$nests[[i]],steps[stepNum]) #Number of foragers to transfer
      moves=whichMoves(nestSet,i=i) #Worst and best cells
      if(moves$move){ #If a move should be made
        #This is not properly taking away foragers from worstNests
        nestSet=moveForagers(nestSet,i,moves) #Moves foragers from cell to cell
      } else if(transfer>1) { #If transfer number is >1 (i.e. not at the end of the list)
        if(verbose) print(with(nestSet$base$nests[[i]],paste0('Finished pass ',stepNum,' of ',
                                length(steps),' for nest ',i,'. Starting pass ',stepNum+1,'.')))
        nestSet$base$nests[[i]]$stepNum=nestSet$base$nests[[i]]$stepNum+1 #Increments step number
        #Creates new best scenario for nest i
        tempBest=makeBest(nestSet$base,whichNest=i,parallel=parallel,cluster=cluster)
        if(!nestSet$base$nests[[i]]$sol) tempWorst=makeWorst(nestSet$base,whichNest=i,parallel=parallel,cluster=cluster) else tempWorst=NA
        nestSet=list(best=tempBest,base=nestSet$base,worst=tempWorst) #Create new scenario set
      } else {done[i]=T} #If transfer is 1, and there's no better deal, distribution has converged
    } #End of FOR loop (nests)
    #If all nests are "done", loop should exit
    nitt=nitt+1 #Increment counter
    if(nitt==iterlim) {
      if(verbose) print('Iteration limit reached')
      break
    }
  } #End of WHILE loop

  #Calculate patch residence times per forager (in seconds)
  for(i in 1:nNests){
    nestSet$base$nests[[i]]$loadingTime=with(nestSet$base,
                                             ifelse(nests[[i]]$n>0,nests[[i]]$L*(nests[[i]]$h+world$S*world$l*nests[[i]]$p_i+world$f)/world$S*world$l,NA))
    nestSet$base$nests[[i]]$travelTime=with(nestSet$base$nests[[i]],ifelse(n>0,(d*(2-beta*(L/L_max)))/(v*(1-beta*(L/L_max))),NA))
    nestSet$base$nests[[i]]$boutLength=with(nestSet$base$nests[[i]],loadingTime+travelTime+H) #Time for 1 complete foraging bout
  }
  if(parallel) {
    stopCluster(cluster) #Stops clusters
    # if(parMethod=='MPI') mpi.finalize() #Cleans MPI states and detaches Rmpi
  }

  if(verbose) print(paste('Simulation ended at ',Sys.time(),
                          '. Elapsed time = ',
                          round(as.numeric(Sys.time()-startTime),3),' ',
                          units(Sys.time()-startTime),
                          '. Final number of iterations = ',
                          nitt,'.',sep=''))
  return(nestSet$base) #Returns world and nests in a list
}
