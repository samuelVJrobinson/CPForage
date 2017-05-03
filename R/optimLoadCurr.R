#Function to calculate currency(-ies) in cell u of nests. Returns NA values if arguments to cells are NA.
optimLoadCurr=function(u,nests,world){
  #Goal: Optimize a vector of L values to produce greatest summed currency in cell u

  #If nest-level arguments are NA or missing, throw an error
  if(any(sapply(nests,function(x) any(lengths(x)==0|is.na(x))))){
    stop('Nest-level arguments are NA or length==0')
  }
  argNames=c("xloc","yloc","n","whatCurr","sol","eps","L_max","v","beta","p_i",
             "h","c_f","c_i","H","d","L","curr")
  if(any(sapply(nests,function(x) any(!argNames %in% names(x))))){
    stop('Nest-level arguments are missing for: ',names(nests)[sapply(nests,function(x) any(!argNames %in% names(x)))])
  }

  emptyNest=sapply(nests,function(x) x$n[u])==0 #Nests in cell U that have no foragers
  #Arguments to feed to optim, which optim feeds to curr
  arglist=list(L_max_i=sapply(nests[!emptyNest],function(x) x$L_max),
               n_i=sapply(nests[!emptyNest],function(x) x$n[u]),
               h_i=sapply(nests[!emptyNest],function(x) x$h[u]),
               p_i=sapply(nests[!emptyNest],function(x) x$p_i),
               f_i=world$f[u],
               d_i=sapply(nests[!emptyNest],function(x) x$d[u]),
               v_i=sapply(nests[!emptyNest],function(x) x$v),
               beta_i=sapply(nests[!emptyNest],function(x) x$beta),
               H_i=sapply(nests[!emptyNest],function(x) x$H),
               c_i=sapply(nests[!emptyNest],function(x) x$c_i),
               c_f=sapply(nests[!emptyNest],function(x) x$c_f),
               whatCurr_i=sapply(nests[!emptyNest],function(x) x$whatCurr),
               mu=world$mu[u],l=world$l[u],e=world$e[u],NumFls=with(world,flDens[u]*cellSize^2),
               patchLev=world$patchLev)

  nestArgs=arglist[c("L_max_i","n_i","h_i","p_i","f_i","d_i","v_i",
                     "beta_i","H_i","c_i","c_f","whatCurr_i")] #Nest-level arguments
  patchArgs=arglist[c('mu','e','NumFls','l')] #Patch-level arguments

  #If anything in the patch-level argument list is NA or <=0, returns NA values - indicates worthless patch
  if(any(sapply(patchArgs,function(x) is.na(x)||x<=0))) {
    #Nest-specific numbers are named after corresponding nest
    return(list('optimL'=setNames(rep(NA,length(nests)),names(nests)),
                #Currency for foraging in empty patches
                'optimCurr'=setNames(sapply(nests,function(x) switch(x$whatCurr,eff=-1,rat=-Inf)),names(nests)),
                'S'=NA))
  }
  #If patch is not worthless, but all nest slots are empty, returns 0 values for L and curr, and 1 for S
  if(sum(emptyNest)==length(emptyNest)){
    return(list('optimL'=setNames(rep(0,length(nests)),names(nests)),
                'optimCurr'=setNames(rep(0,length(nests)),names(nests)),
                'S'=1)) #No competition in completely empty cells
  }


  startL=sapply(nests[!emptyNest],function(x) x$L[u]) #Starting vector for L-values
  startL[is.na(startL)]=0 #If there are any nests with NA L-values, sets L-value to zero (arglist already checked for NAs, so this means that the cell hasn't been used yet)
  startL[startL!=0]=0 #TEMPORARY: ALWAYS USES 0 AS STARTING VALUE FOR LOAD
  #L values that produce highest sum of currency
  optimL=do.call(optim,c(list(par=startL,fn=curr,method='L-BFGS-B',lower=rep(0,length(startL)),
                              upper=sapply(nests[!emptyNest],function(x) x$L_max),control=list(fnscale=-1)),arglist))$par
  #Currency for each nest, and S-value for the cell
  currencyS=do.call(curr,c(list(L=optimL,sumAll=F),arglist)) #Named vector of currency and S-values
  optimCurr=currencyS[1:length(currencyS)-1]
  S=currencyS[length(currencyS)]
  if(sum(emptyNest)>0) { #If some nests had no foragers in cell U
    tempL=tempCurr=setNames(rep(0,length(nests)),names(nests)) #Vector of zeros
    tempL[names(nests)==names(optimL)]=optimL #Adds calculated values alongside zeros
    tempCurr[names(nests)==names(optimCurr)]=optimCurr
    optimL=tempL #Overwrites old vector
    optimCurr=tempCurr
  }
  # Return all results together in one list
  resultList=list('optimL'=optimL,'optimCurr'=optimCurr,'S'=S)
  return(resultList)
}
