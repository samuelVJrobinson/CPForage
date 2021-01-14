#Convert CPForage model to dataframe
#' Convert CPForage model (list structure) to dataframe.
#'
#' @param modelRun CPForage model
#'
#' @return Dataframe containing spatial data from CPForage simulation
#'
#' @examples
#'
#' #Create test world for run
#' nu_i<-0.15/3600 #Nectar production/hr for a single flower
#' flDens<-200 #Flower density/m2 - represents early canola bloom
#' e_i<-14.35 #Energetic value/unit
#' l_i<-0.87 #Canola standing crop (0.87uL)
#' f_i<-3 #Inter-flower flight time
#' #World structure
#' cellSize<-10 #10m cells (100m^2)
#' worldSize<-120 #120x120m field (100x100m field with 10m buffer zone with weeds in it)
#' world1<-list(mu=matrix(1/3600,worldSize/cellSize,worldSize/cellSize),  #Weeds produce 1uL/hr
#'              flDens=matrix(3*(cellSize^2),worldSize/cellSize,worldSize/cellSize), #1 fls/m2
#'              e=matrix(8,worldSize/cellSize,worldSize/cellSize), #worth about 75% of canola
#'              l=matrix(3,worldSize/cellSize,worldSize/cellSize), #max nectar = 3uL
#'              f=matrix(5,worldSize/cellSize,worldSize/cellSize), #5 second flight time b/w fls
#'              alphaVal=matrix(0.013,worldSize/cellSize,worldSize/cellSize), #Alpha value (J/s*uL)
#'              cellSize=cellSize,
#'              forageType='random') #Competition for flowers within patch
#' world1$mu[c(2:11),c(2:11)]<-nu_i #Per-flower nectar production in
#' # canola-filled cells
#' world1$flDens[c(2:11),c(2:11)]<-flDens*cellSize^2 #Flower number per cell
#' world1$e[c(2:11),c(2:11)]<-e_i #Energy production in canola-filled cells
#' world1$l[c(2:11),c(2:11)]<-l_i #Standing crop in cells with no competition
#' world1$f[c(2:11),c(2:11)]<-f_i #Inter-flower flight time world1$patchLev
#'
#' #Constants for foragers
#' honeybeeConstants<-list(L_max=59.5, #Max load capacity (uL)
#'                         v=7.8, #Velocity (m/s) - Unloaded flight speed
#'                         betaVal=0.102/59.5, #Proportion reduction in completely loaded flight speed (1-v/v_l)
#'                         p_i=1, # Max loading rate (uL/s)
#'                         h=1.5, #Handling time per flower (s)
#'                         c_f=0.05, #Unloaded flight energetic cost (J/s)
#'                         c_i=0.0042, #Cost of non-flying activity
#'                         H=100) #Time spent in the hive (s)
#'
#' #Nest structure (social rate maximizers)
#' nests1<-list(xloc=1,yloc=1,n=1000,whatCurr='rat',sol=FALSE,
#'              constants=honeybeeConstants,eps=0,steps=c(50,5,1))
#'
#' #Run model
#' testOutput1<-forageMod(world1,nests1,2000,verbose=FALSE,parallel=FALSE)
#'
#' nests_df <- cpf2df(nests1) #Convert to dataframe
#' head(nests_df)
#'

cpf2df <- function(modelRun){

  xdim <- ncol(modelRun$nests$n) #Dimensions of world matrix
  ydim <- nrow(modelRun$nests$n)

  d <- data.frame(x=rep(1:xdim,each=ydim),y=rep(1:ydim,times=xdim)) #Cell locations

  #Terms from world
  worldTerms <- c('mu','flDens','e','l','f','alphaVal','S')
  d <- cbind(d,sapply(worldTerms,function(x) as.vector(modelRun$world[[x]])))

  #Terms from nest
  nestTerms <- c('n','L','curr','d','travelTime','loadingTime','boutLength')
  d <- cbind(d,sapply(nestTerms,function(x) as.vector(modelRun$nests[[x]])))

  return(d)
}
