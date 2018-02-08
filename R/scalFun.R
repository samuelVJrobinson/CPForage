#' scalFun
#'
#' Function to find  S (scaling term), given Load, Forager number, Max intake rate, travel time, and nectar production
#'
#'@param mu Per-flower nectar production (\eqn{\muL}/s)
#'@param l Maximum standing crop per flower (\eqn{\muL})
#'@param NumFls Number of flowers per patch
#'@param L_i Load size (\eqn{\muL})
#'@param L_max_i Maximum load size (\eqn{\muL})#'
#'@param n_i Number of foragers in the cell
#'@param h_i Handling time per flower (s)
#'@param p_i Licking speed for nectar (\eqn{\muL/s})
#'@param f_i Flight time between flowers (s)
#'@param d_i Distance from hive (m)
#'@param v_i Unloaded flight speed from hive (m/s)
#'@param beta_i Decrease in flight speed with load (m/s\eqn{\muL})
#'@param H_i Time spent inside hive (s)
#'@param forageType Type of foraging within patch (see details)
#'
#' @return S-value for patch
#'
#' @details \code{forageType} can be one of the following:
#' \itemize{
#' \item \code{random}: random foraging (Possingham 1988)
#' \item \code{multiNN}: nearest-neighbour foraging
#' \item \code{omniscient}: flowers are shared optimally between foragers. Non-saturating nectar production.
#' \item \code{random_nonsat}: random foraging with non-saturating nectar production (ceiling at \code{l}).
# \item \code{singleNN_nonsat}: nearest-neighbour foraging with non-saturating nectar production, assuming all visits done by one forager
#'
#' }
#' @examples
#' #Test
#' mu=0.3/3600
#' l=1
#' NumFls=520 #520 flowers
#' L_i=59.5
#' L_max_i=59.5
#' n_i=10
#' h_i=1.5
#' p_i=1
#' f_i=0.86
#' d_i=100
#' v_i=7.8
#' beta_i=0.102
#' H_i=100
#' forageType='omniscient'
#'
#' scalFun(mu,l,NumFls,L_i,L_max_i,n_i,h_i,p_i,f_i,
#'                 d_i,v_i,beta_i,H_i,forageType)
scalFun=function(mu=NULL,l=NULL,NumFls=NULL,L_i=NULL,L_max_i=NULL,n_i=NULL,h_i=NULL,p_i=NULL,f_i=NULL,d_i=NULL,v_i=NULL,beta_i=NULL,H_i=NULL,forageType=NULL){

  #This function is intended to be run for a single cell (scalar arguments for
  #patch-level properties: mu,l,NumFls), but can return S-values for multiple
  #nexts (scalar or vector arguments for all other properties)

  #Argument checking
  singleArgs=list(mu=mu,l=l,NumFls=NumFls) #Patch arguments (mu,l,NumFls)
  #If any arguments are missing, throw an error
  if(any(is.na(singleArgs))|any(sapply(singleArgs,is.null))){
    stop('Patch-level arguments ',names(singleArgs)[sapply(singleArgs, function(x) is.null(x)||is.na(x))],' missing')
  }
  singArgLen=sapply(singleArgs,length) #Length of patch-level arguments
  if(sum(singArgLen>1)>0){ #If patch arguments are longer than 1
    stop(names(singArgLen)[singArgLen>1],'longer than 1. Patch arguments must be a scalar.')
  }
  #If nectar production is 0, competition isn't defined.
  if(singleArgs$mu==0|singleArgs$l==0|singleArgs$NumFls==0) return(NA)

  #If L_i is less than 0 (negative competition), S = 1 (optimizer occasionally has to use
  #slightly negative L values)
  if(L_i<0) return(1)

  #Forager (nest-level) arguments (Load,MaxLoad,number,Handling Time,Licking
  #Rate,b/w flower flight time,distance,beta,hive time)
  Args=list(L_i=L_i,L_max_i=L_max_i,n_i=n_i,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i)
  argLen=sapply(Args,length) #Length of forager arguments
  argSign=sapply(c(singleArgs,Args),function(x) sum(x<0))>0 #Are any arguments less than 0?
  argSign=argSign[names(argSign)!='L_i'] #Removes L_i from argSign
  #If any arguments (other than L_i) are <0
  if(any(argSign)){
    stop(names(argSign)[argSign],' < 0')
  } else if(any(argLen==0)|any(singArgLen==0)) {#If there are any zero-length arguments (i.e. not provided)
    stop(paste(names(singArgLen)[singArgLen==0],names(argLen)[argLen==0],'arguments not provided'))
  } else if(any(argLen>1)) stop("Too many nest arguments. CPforage can't handle scenarios with >1 nest.")

  #For omniscient and random foragers there is a definite solution for S|L,other
  #parameters. Nearest neighbour foraging is nonlinear, so solution must be
  #found iteratively

  S=switch(forageType,
           #Omniscient
          omniscient={
            X = -(L_i*NumFls*(h_i+f_i)*mu*v_i*(L_i*beta_i-L_max_i))/(l*(L_i^2*NumFls*mu*p_i*v_i*beta_i-
            L_i^2*n_i*v_i*beta_i+H_i*L_i*NumFls*mu*v_i*beta_i+L_i*NumFls*d_i*mu*beta_i-L_i*L_max_i*NumFls*mu*
            p_i*v_i+L_i*L_max_i*n_i*v_i-H_i*L_max_i*NumFls*mu*v_i-2*L_max_i*NumFls*d_i*mu))
            ifelse(X<=0,1,X) #Sets S to 1 if less than 0 (can't have negative competition)
            },
           #Random - Possingham 1988: 1/(D_lamda*l+1)=S
          random=-(sqrt(((mu^2*L_i^4*NumFls^2*l^2*p_i^2+(-2*mu*L_i^4*NumFls*l^2*n_i+2*mu^2*H_i*L_i^3*NumFls^2*
            l^2+(2*mu^2*L_i^4*NumFls^2*h_i+2*mu^2*L_i^4*NumFls^2*f_i)*l)*p_i+L_i^4*l^2*n_i^2+((2*mu*L_i^4*
            NumFls*h_i+2*mu*L_i^4*NumFls*f_i)*l-2*mu*H_i*L_i^3*NumFls*l^2)*n_i+mu^2*H_i^2*L_i^2*NumFls^2*l^2+
            (2*mu^2*H_i*L_i^3*NumFls^2*h_i+2*mu^2*H_i*L_i^3*NumFls^2*f_i)*l+mu^2*L_i^4*NumFls^2*h_i^2+2*mu^2*
            L_i^4*NumFls^2*f_i*h_i+mu^2*L_i^4*NumFls^2*f_i^2)*v_i^2+(2*mu^2*L_i^3*NumFls^2*d_i*l^2*p_i-2*mu*
            L_i^3*NumFls*d_i*l^2*n_i+2*mu^2*H_i*L_i^2*NumFls^2*d_i*l^2+(2*mu^2*L_i^3*NumFls^2*d_i*h_i+2*mu^2*
            L_i^3*NumFls^2*d_i*f_i)*l)*v_i+mu^2*L_i^2*NumFls^2*d_i^2*l^2)*beta_i^2+((-2*mu^2*L_i^3*L_max_i*NumFls^2*
            l^2*p_i^2+(4*mu*L_i^3*L_max_i*NumFls*l^2*n_i-4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*l^2+(-4*mu^2*L_i^3*
            L_max_i*NumFls^2*h_i-4*mu^2*L_i^3*L_max_i*NumFls^2*f_i)*l)*p_i-2*L_i^3*L_max_i*l^2*n_i^2+(4*mu*H_i*
            L_i^2*L_max_i*NumFls*l^2+(-4*mu*L_i^3*L_max_i*NumFls*h_i-4*mu*L_i^3*L_max_i*NumFls*f_i)*l)*n_i-2*
            mu^2*H_i^2*L_i*L_max_i*NumFls^2*l^2+(-4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*h_i-4*mu^2*H_i*L_i^2*L_max_i*
            NumFls^2*f_i)*l-2*mu^2*L_i^3*L_max_i*NumFls^2*h_i^2-4*mu^2*L_i^3*L_max_i*NumFls^2*f_i*h_i-2*mu^2*
            L_i^3*L_max_i*NumFls^2*f_i^2)*v_i^2+(-6*mu^2*L_i^2*L_max_i*NumFls^2*d_i*l^2*p_i+6*mu*L_i^2*L_max_i*
            NumFls*d_i*l^2*n_i-6*mu^2*H_i*L_i*L_max_i*NumFls^2*d_i*l^2+(-6*mu^2*L_i^2*L_max_i*NumFls^2*d_i*h_i-6*
            mu^2*L_i^2*L_max_i*NumFls^2*d_i*f_i)*l)*v_i-4*mu^2*L_i*L_max_i*NumFls^2*d_i^2*l^2)*beta_i+(mu^2*L_i^2*
            L_max_i^2*NumFls^2*l^2*p_i^2+(-2*mu*L_i^2*L_max_i^2*NumFls*l^2*n_i+2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*
            l^2+(2*mu^2*L_i^2*L_max_i^2*NumFls^2*h_i+2*mu^2*L_i^2*L_max_i^2*NumFls^2*f_i)*l)*p_i+L_i^2*L_max_i^2*
            l^2*n_i^2+((2*mu*L_i^2*L_max_i^2*NumFls*h_i+2*mu*L_i^2*L_max_i^2*NumFls*f_i)*l-2*mu*H_i*L_i*L_max_i^2*
            NumFls*l^2)*n_i+mu^2*H_i^2*L_max_i^2*NumFls^2*l^2+(2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*h_i+2*mu^2*H_i*
            L_i*L_max_i^2*NumFls^2*f_i)*l+mu^2**NumFls^2*h_i^2+2*mu^2*L_i^2*L_max_i^2*NumFls^2*f_i*h_i+mu^2*
            L_i^2*L_max_i^2*NumFls^2*f_i^2)*v_i^2+(4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*l^2*p_i-4*mu*L_i*L_max_i^2*
            NumFls*d_i*l^2*n_i+4*mu^2*H_i*L_max_i^2*NumFls^2*d_i*l^2+(4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*h_i+4*mu^2*
            L_i*L_max_i^2*NumFls^2*d_i*f_i)*l)*v_i+4*mu^2*L_max_i^2*NumFls^2*d_i^2*l^2)+((-mu*L_i^2*NumFls*l*p_i+
            L_i^2*l*n_i-mu*H_i*L_i*NumFls*l+mu*L_i^2*NumFls*h_i+mu*L_i^2*NumFls*f_i)*v_i-mu*L_i*NumFls*d_i*l)*
            beta_i+(mu*L_i*L_max_i*NumFls*l*p_i-L_i*L_max_i*l*n_i+mu*H_i*L_max_i*NumFls*l-mu*L_i*L_max_i*NumFls*
            h_i-mu*L_i*L_max_i*NumFls*f_i)*v_i+2*mu*L_max_i*NumFls*d_i*l)/(((2*mu*L_i^2*NumFls*l*p_i+2*mu*H_i*
            L_i*NumFls*l)*v_i+2*mu*L_i*NumFls*d_i*l)*beta_i+(-2*mu*L_i*L_max_i*NumFls*l*p_i-2*mu*H_i*L_max_i*
            NumFls*l)*v_i-4*mu*L_max_i*NumFls*d_i*l),
           #Random foraging, with non-saturating nectar production
          random_nonsat=do.call(optimize,list(f=nonlinearS,interval=c(0,2),
                A=0.6100842,B=1.2001725,C=0,
                mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=n_i,h_i=h_i,p_i=p_i,
                f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i))$minimum,
          ## Single nearest-neighbour foraging, with non-saturating nectar production -not that useful
          # singleNN_nonsat=do.call(optimize,list(f=nonlinearS,interval=c(0,2),
          #       A=5.373864e-01 + 1.968079e+04*(mu/l) - 1.060072e+08*(mu/l)^2,
          #       B=1.546423e+00 - 1.427107e+03*(mu/l) - 2.414069e+07*(mu/l)^2,
          #       C=6.870373e-03 + 3.292027e+03*(mu/l) - 1.173807e+07*(mu/l)^2,
          #       mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=n_i,h_i=h_i,p_i=p_i,
          #       f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i))$minimum,
          multiNN=do.call(optimize,list(f=nonlinearS_2,interval=c(0,1),
                a0=0.900152061, #Taken from simulation results
                a1=0.005108424,
                b1=0.001209028,
                mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=n_i,h_i=h_i,p_i=p_i,
                f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i))$minimum,
          stop(paste0('Foraging type: "',forageType,'" not found.'))
          )

  #If S is greater than 1, sets it back to 1 (flowers can't give more than their max).
  if(S>1) return(1) else return(S)
}
