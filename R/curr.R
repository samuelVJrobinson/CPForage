#'Extract currency (or currencies) from a cell.
#'
#Function to calculate currency(-ies) in a given cell. Called by \code{optimLoadCurr}
#'
#'@param L_i Load size (\eqn{\muL})
#'@param L_max_i Maximum load size (\eqn{\muL})
#'@param n_i Number of foragers in the cell
#'@param h_i Handling time per flower (s)
#'@param p_i Licking speed for nectar (\eqn{\muL/s})
#'@param f_i Flight time between flowers (s)
#'@param d_i Distance from hive (m)
#'@param v_i Unloaded flight speed from hive (m/s)
#'@param beta_i Decrease in flight speed with load (m/s\eqn{\muL})
#'@param H_i Time spent inside hive (s)
#'@param c_i Cost of non-flying behaviour (J/s)
#'@param c_f Cost of flight (J/s)
#'@param whatCurr_i Currency to use. Must be either "rat" (net rate) or "eff" (efficiency)
#'@param mu Per-flower nectar production (\eqn{\muL}/s)
#'@param l Maximum standing crop per flower (\eqn{\muL})
#'@param e Energetic value of nectar (J/\eqn{\muL})
#'@param NumFls Number of flowers per patch
#'@param sumAll Should currencies for the cell be summed? (useful for optimization of L) If not, returns a vector of currencies.

#'@return Summed currency, or vector of currencies
#'
#'@examples
#'#Parameters
#'params=list(L_i=59.5,
#'L_max_i=59.5,
#'n_i=10,
#'h_i=1.5,
#'p_i=1,
#'f_i=0.86,
#'d_i=100,
#'v_i=7.8,
#'beta_i=0.102,
#'H_i=100,
#'c_f=0.05,
#'c_i=0.0042,
#'mu=0.3/3600,
#'l=1,
#'e=14,
#'NumFls=520*(10^2))
#'
#'#Currency for 1 nest
#'with(params,
#'curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#'c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,patchLev=F)
#')
#'
#'#Summed currency for 2 nests
#'with(params,
#'curr(L_i,L_max_i,n_i,h_i,p_i,f_i,c(d_i,300),v_i,beta_i,H_i,
#'c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,patchLev=F)
#')
#'
#'#Individual currencies for 2 nests
#'with(params,
#'curr(L_i,L_max_i,n_i,h_i,p_i,f_i,c(d_i,300),v_i,beta_i,H_i,
#'c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,patchLev=F)
#')

curr=function(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,c_i,c_f,whatCurr_i,mu,l,e,NumFls,sumAll=T,patchLev=F){
  Stemp=scalFun(mu=mu, l=l, L_i=L_i, L_max=L_max_i, n_i=n_i, #Scaling function for cell
                h_i=h_i, p_i=p_i, f_i=f_i, d_i=d_i, v_i=v_i, beta_i=beta_i, H_i=H_i,NumFls=NumFls,patchLev=patchLev)
  #Vector of currencies for cell, using getCurr to return appropriate nest-specific currency
  currs=mapply(getCurr,whatCurr=whatCurr_i, L=L_i,L_max=L_max_i, e=e, d=d_i, v=v_i, h=h_i,
               f=f_i, l=l, p_i=p_i, c_i=c_i, c_f=c_f, H=H_i, beta=beta_i,Stemp)
  #If sumAll is T, add all currencies together and return (useful for optimization of L).
  #Else, returns vector of currencies and S value for that cell (useful for generally returning currency)
  if(sumAll) return (sum(currs)) else return(setNames(c(currs,Stemp),c(names(currs),'S')))
}
