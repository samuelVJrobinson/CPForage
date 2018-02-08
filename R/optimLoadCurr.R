#' Optimal currency and load
#'
#' Function to calculate optimal currency and load in cell u. Returns NA values
#' if arguments to cells are NA. This is the "workhorse" function of CPForage
#'
#' @param u Address of cells in Scenario to calculate
#' @param scenario Scenario to use
#'
#' @return List of optimal currency and load values for each cell in \code{u}
#'
#' @examples
#'

optimLoadCurr=function(u,scenario){
  #Goal: Optimize a vector of L values to produce greatest summed currency in cell u
  nests=scenario$nests #Unpacks scenario
  world=scenario$world
  #If nest-level arguments are NA or missing, throw an error
  if(any(sapply(nests,function(x) any(lengths(x)==0|is.na(x))))){
    stop('Nest-level arguments are NA or length==0')
  }
  argNames=c("xloc","yloc","n","whatCurr","sol","eps","L_max","v","beta","p_i",
             "h","c_f","c_i","H","d","L","curr")
  if(any(sapply(nests,function(x) any(!argNames %in% names(x))))){
    stop('Nest-level arguments are missing for: ',paste(names(nests)[sapply(nests,function(x)
      any(!argNames %in% names(x)))]))
  }
  #If there are too many cells provided
  if(length(u)>1) {
    stop('Too many cells provided. Use lapply to pass cells. e.g. lapply(use,optimLoadCurr,scenario=scenario)')
  }

  emptyNest=nests[[1]]$n[u]<1 #Is cell u empty?

  #If cell u is empty, returns NA values for L, 0 for curr, and 1 for S
  if(emptyNest){
    return(list('optimL'=setNames(rep(NA,length(nests)),names(nests)),
                'optimCurr'=setNames(rep(0,length(nests)),names(nests)),
                'S'=1)) #No competition in completely empty cells
  }

  #Arguments to feed to optim, which optim feeds to curr
  arglist=list(L_max_i=nests[[1]]$L_max,n_i=nests[[1]]$n[u],
               h_i=nests[[1]]$h[u],
               p_i=nests[[1]]$p_i,
               f_i=world$f[u],
               d_i=nests[[1]]$d[u],
               v_i=nests[[1]]$v,
               beta_i=nests[[1]]$beta,
               H_i=nests[[1]]$H,
               c_i=nests[[1]]$c_i,
               c_f=nests[[1]]$c_f,
               whatCurr_i=nests[[1]]$whatCurr,
               mu=world$mu[u],l=world$l[u],e=world$e[u],NumFls=with(world,flDens[u]*cellSize^2),
               forageType=world$forageType)


  #Nest-level arguments (one for each nest involved)
  nestArgs=arglist[c("L_max_i","n_i","p_i","f_i","d_i","v_i",
                     "beta_i","H_i","c_i","c_f","whatCurr_i","forageType")]
  #Patch-level arguments (only one for the patch)
  patchArgs=arglist[c('mu','e','NumFls','l','h_i')]

  #Are any nest-level arguments NA or nonexistant?
  if(any(sapply(nestArgs,function(x) any(lengths(x)==0|is.na(x))))){
    stop('Nest-level arguments ',paste(names(nestArgs)[sapply(nestArgs,function(x) any(lengths(x)==0|is.na(x)))]),' are NA or length==0. Are all dimensions equal?')
  }
  #Are any patch-level arguments nonexistent?
  if(any(sapply(patchArgs,function(x) any(lengths(x)==0)))) {
    stop('Patch-level arguments ',paste(names(patchArgs)[sapply(patchArgs,function(x) any(lengths(x)==0|is.na(x)))]),' are missing (length==0). Are all dimensions equal?')
  }

  #If anything in the patch-level argument list is NA or <=0, returns NA values - indicates worthless patch
  if(any(sapply(patchArgs,function(x) is.na(x)||x<=0))) {
    #Nest-specific numbers are named after corresponding nest
    return(list('optimL'=setNames(rep(NA,length(nests)),names(nests)),
                #Currency for foraging in empty patches
                'optimCurr'=setNames(sapply(nests,function(x) switch(x$whatCurr,eff=-1,rat=-Inf)),names(nests)),
                'S'=NA))
  }
  startL=0 #ALWAYS USES 0 AS STARTING VALUE FOR LOAD
  #L value and maximized currency value
  optimL=do.call(optimize,c(list(f=curr,interval=c(0,nests[[1]]$L_max),maximum=T),arglist))

  # #Multi-nest version
  # optimL=do.call(optim,c(list(par=startL,fn=curr,method='L-BFGS-B',lower=0,upper=nests[[1]]$L_max,
  #                             control=list(fnscale=-1)),arglist))$par

  #Best currency given optimum load, and S-value for the cell
  currencyS=do.call(curr,c(list(L=optimL$maximum,sumAll=F),arglist)) #Named vector of currency and S-values
  optimCurr=currencyS[[1]]
  S=currencyS[[2]]

  # Return all results together in one list
  resultList=list('optimL'=optimL$maximum,'optimCurr'=optimCurr,'S'=S)
  return(resultList)
}
