#' Optimal currency and load
#'
#' Function to calculate optimal currency and load in cell u. Returns NA values
#' if arguments to cells are NA. This is the "workhorse" function of CPForage.
#'
#' @param u Address of cells in Scenario to calculate
#' @param scenario Scenario to use
#'
#' @return List of optimal currency and load values for each cell in \code{u}
#'
#' @examples
#'

optimLoadCurr <- function(u,scenario){
  #Goal: Optimize a vector of L values to produce greatest summed currency in cell u
  nests <- scenario$nests #Unpacks scenario
  world <- scenario$world

  #Argument checking:
  #If nest-level arguments are NA or missing, throw an error
  if(any(lengths(nests)==0|is.na(nests))){
    stop('Nest-level arguments are NA or length==0')
  }
  #If any arguments are missing from nests
  argNames <- c("xloc","yloc","n","whatCurr","sol","eps","L_max","v","beta","p_i",
             "h","c_f","c_i","H","d","L","curr")
  if(any(!argNames %in% names(nests))){
    stop('Nest-level arguments are missing for: ',paste(names(nests)[sapply(nests,function(x)
      any(!argNames %in% names(x)))],sep=','))
  }
  #If there are too many cells provided
  if(length(u)>1) {
    stop('Too many cells provided. Use lapply to pass cells. e.g. lapply(use,optimLoadCurr,scenario=scenario)')
  }

  #If cell u has no foragers, returns NA values for L, 0 for curr, and 1 for S
  if(nests$n[u]<1){
    return(list('optimL'=NA,'optimCurr'=0,'S'=1)) #No competition in unoccupied cells
  }

    #Arguments to feed to optim, which optim feeds to curr
  arglist=list(L_max_i=nests$L_max,n_i=nests$n[u],
               h_i=nests$h[u],
               p_i=nests$p_i,
               d_i=nests$d[u],
               v_i=nests$v,
               beta_i=nests$beta,
               H_i=nests$H,
               c_i=nests$c_i,
               c_f=nests$c_f,
               whatCurr_i=nests$whatCurr,
               mu=world$mu[u],l=world$l[u],e=world$e[u],NumFls=world$flDens[u],
               f_i=world$f[u],forageType=world$forageType,
               alphaVal=world$alphaVal[u])

  #Nest-level arguments (one for each nest involved)
  nestArgs <- arglist[c("L_max_i","n_i","p_i","f_i","d_i","v_i",
                     "beta_i","H_i","c_i","c_f","whatCurr_i","forageType")]
  #Patch-level arguments (only one for the patch)
  patchArgs <- arglist[c('mu','e','NumFls','l','h_i','alphaVal')]

  #Are any nest-level arguments NA or nonexistant?
  if(any(lengths(nestArgs)==0|is.na(nestArgs))){
    stop('Nest-level arguments ',paste(names(nestArgs)[any(lengths(nestArgs)==0|is.na(nestArgs))]),
         ' are NA or length==0. Are all dimensions equal?')
  }
  #Are any patch-level arguments nonexistent?
  if(any(lengths(patchArgs)==0)) {
    stop('Patch-level arguments ',paste(names(patchArgs)[any(lengths(patchArgs)==0|is.na(patchArgs))]),
         ' are missing (length==0). Are all dimensions equal?')
  }

  #If anything in the patch-level argument list is NA or <=0, returns NA values - indicates worthless patch
  if(any(is.na(patchArgs)|patchArgs<=0)) {
    return(list('optimL'=NA,'optimCurr'= switch(nests$whatCurr,eff=-1,rat=-Inf),'S'=NA))
  }

  startL <- 0 #Use zero as the starting value for load

  #L value and maximized currency value - tolerance needs to be < ~1e-7
  optimL <- do.call(optimize,c(list(f=curr,interval=c(0,nests$L_max),maximum=TRUE,tol=1e-10),arglist))

  #Best currency given optimum load, and S-value for the cell
  #NOTE: this works for both solitary and social, because it calculates (currency | n); n is dealt with elsewhere
  currencyS <- do.call(curr,c(list(L=optimL$maximum,sumAll=FALSE),arglist)) #Named vector of currency and S-values
  optimCurr <- currencyS[[1]]
  S <- currencyS[[2]]

  # Return all results together in one list
  resultList <- list('optimL'=optimL$maximum,'optimCurr'=optimCurr,'S'=S)
  return(resultList)
}
