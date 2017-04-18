#Function to calculate currency(-ies) in cell u of nests. Returns NA values if arguments to cells are NA.
optimLoadCurr=function(u,nests,world){
  #Optimize a vector of L values to produce greatest summed currency in cell u
  emptyNest=sapply(nests,function(x) x$n[u])==0 #Nests in cell U that have no foragers
  if(length(emptyNest)==0){ #If all nests are empty, returns 0 values for L and curr, and 1 for S
    return(list('optimL'=setNames(rep(0,length(nests)),names(nests)),
                'optimCurr'=setNames(rep(0,length(nests)),names(nests)),
                'S'=1)) #No competition in completely empty cells
  }
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
               mu=world$mu[u],l=world$l[u],e=world$e[u],NumFls=with(world,flDens*cellSize^2))
  if(sum(sapply(arglist,is.na))>0) { #If anything in the argument list is NA, returns NA values
    #Nest-specific numbers are named after corresponding nest
    return(list('optimL'=setNames(rep(NA,length(nests)),names(nests)),
                #Currency for foraging in empty patches
                'optimCurr'=setNames(sapply(nests,function(x) switch(x$whatCurr,eff=-1,rat=-Inf)),names(nests)),
                'S'=NA))
  }
  startL=sapply(nests[!emptyNest],function(x) x$L[u]) #Starting vector for L-values
  startL[is.na(startL)]=0 #If there are any nests with NA L-values, sets L-value to zero (arglist already checked for NAs, so this means that the cell hasn't been used yet)
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
