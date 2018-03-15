#'Extract currency (or currencies) from a cell.
#'
#'Function to calculate currency(-ies) in a given cell. Called by \code{optimLoadCurr}
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
#'@param forageType Foraging style (see details below).
#'
#'@details \code{forageType} can be one of the following:
#' \itemize{
#' \item \code{random}: random foraging (Possingham 1988)
#' \item \code{nn}: nearest-neighbour foraging
#' \item \code{omniscient}: flowers are shared optimally between foragers. Non-saturating nectar production.
#' \item \code{random_nonsat}: random foraging with non-saturating nectar production (ceiling at \code{l}).
#' }

#'@return Summed currency, or vector of currencies
#'
#'@examples
#'#Parameters
#'params=list(L_i=59.5,
#'            L_max_i=59.5,
#'            n_i=10,
#'            h_i=1.5,
#'            p_i=1,
#'            f_i=0.86,
#'            d_i=100,
#'            v_i=7.8,
#'            beta_i=0.102,
#'            H_i=100,
#'            c_f=0.05,
#'            c_i=0.0042,
#'            mu=0.3/3600,
#'            l=1,
#'            e=14.35,
#'            NumFls=520*(10^2))
#'
#'temp1 <- with(params,
#'     curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#'              c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,
#'              forageType='omniscient'))
#'
#'temp1[['eff']] #Efficiency within cell
#'temp1[['S']] #Reduction in per-flower value


curr=function(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,c_i,c_f,whatCurr_i,mu,l,e,NumFls,sumAll=T,forageType='random'){
  #Calculate S (competition term)
  if(L_i<=0){ #If Load is less than zero, S=1 (impossible to give back to patch)
    S=1
  } else {
    S=switch(forageType,
             #Omniscient - all flowers used optimally
             omniscient={
               X = -(L_i*NumFls*(h_i+f_i)*mu*v_i*(L_i*beta_i-L_max_i))/(l*(L_i^2*NumFls*mu*p_i*v_i*beta_i-
               L_i^2*n_i*v_i*beta_i+H_i*L_i*NumFls*mu*v_i*beta_i+L_i*NumFls*d_i*mu*beta_i-L_i*L_max_i*NumFls*mu*
               p_i*v_i+L_i*L_max_i*n_i*v_i-H_i*L_max_i*NumFls*mu*v_i-2*L_max_i*NumFls*d_i*mu))
               ifelse(X<0,1,X) #Sets S to 1 if less than 0 (can't have negative competition)
               },
             #Random - Possingham 1988: 1/(D_lamda*l+1)=S
             random={-(sqrt(((mu^2*L_i^4*NumFls^2*l^2*p_i^2+(-2*mu*L_i^4*NumFls*l^2*n_i+2*mu^2*H_i*L_i^3*NumFls^2*
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
               NumFls*l)*v_i-4*mu*L_max_i*NumFls*d_i*l)},
             #Random foraging, with non-saturating nectar production
             random_nonsat={do.call(optimize,list(f=nonlinearS,interval=c(0,2),
               A=0.6100842,B=1.2001725,C=0,
               mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=n_i,h_i=h_i,p_i=p_i,
               f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i))$minimum},
             #Nearest-neighbour foraging. Used same as random, as simulation results were very similar.
             nn={Sval=-(sqrt(((mu^2*L_i^4*NumFls^2*l^2*p_i^2+(-2*mu*L_i^4*NumFls*l^2*n_i+2*mu^2*H_i*L_i^3*NumFls^2*
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
                 NumFls*l)*v_i-4*mu*L_max_i*NumFls*d_i*l)
                 DlambdaS=(L_i*n_i)/ #Scaled Dlambda value
                  (mu*NumFls*Sval*((d_i*(2*L_max_i-L_i*beta_i))/(v_i*(L_max_i-L_i*beta_i))+(L_i*(Sval*l*p_i+h_i+f_i))/(Sval*l)+H_i))
                 #If there are less than 5 foragers operating in the patch and DlambdaS<0.5, uses different set of behaviour
                 if(n_i<5 & DlambdaS<0.6){
                   slope=-0.76906
                   brPoint=0.07669
                   Sval=-(sqrt((((mu^2*L_i^4*NumFls^2*brPoint^2*l^2*p_i^2+2*mu^2*H_i*L_i^3*NumFls^2*brPoint^2*l^2*p_i+mu^2*
                    H_i^2*L_i^2*NumFls^2*brPoint^2*l^2)*slope^2+(-2*mu^2*L_i^4*NumFls^2*brPoint*l^2*p_i^2+(4*mu*L_i^4*NumFls*
                    l^2*n_i-4*mu^2*H_i*L_i^3*NumFls^2*brPoint*l^2+(-2*mu^2*L_i^4*NumFls^2*brPoint*h_i-2*mu^2*L_i^4*NumFls^2*
                    brPoint*f_i)*l)*p_i+4*mu*H_i*L_i^3*NumFls*l^2*n_i-2*mu^2*H_i^2*L_i^2*NumFls^2*brPoint*l^2+(-2*mu^2*H_i*L_i^3*
                    NumFls^2*brPoint*h_i-2*mu^2*H_i*L_i^3*NumFls^2*brPoint*f_i)*l)*slope+mu^2*L_i^4*NumFls^2*l^2*p_i^2+(2*mu^2*H_i*
                    L_i^3*NumFls^2*l^2+(2*mu^2*L_i^4*NumFls^2*h_i+2*mu^2*L_i^4*NumFls^2*f_i)*l)*p_i+mu^2*H_i^2*L_i^2*NumFls^2*l^2+
                    (2*mu^2*H_i*L_i^3*NumFls^2*h_i+2*mu^2*H_i*L_i^3*NumFls^2*f_i)*l+mu^2*L_i^4*NumFls^2*h_i^2+2*mu^2*L_i^4*
                    NumFls^2*f_i*h_i+mu^2*L_i^4*NumFls^2*f_i^2)*v_i^2+((2*mu^2*L_i^3*NumFls^2*brPoint^2*d_i*l^2*p_i+2*mu^2*
                    H_i*L_i^2*NumFls^2*brPoint^2*d_i*l^2)*slope^2+(-4*mu^2*L_i^3*NumFls^2*brPoint*d_i*l^2*p_i+4*mu*L_i^3*NumFls*
                    d_i*l^2*n_i-4*mu^2*H_i*L_i^2*NumFls^2*brPoint*d_i*l^2+(-2*mu^2*L_i^3*NumFls^2*brPoint*d_i*h_i-2*mu^2*
                    L_i^3*NumFls^2*brPoint*d_i*f_i)*l)*slope+2*mu^2*L_i^3*NumFls^2*d_i*l^2*p_i+2*mu^2*H_i*L_i^2*NumFls^2*d_i*
                    l^2+(2*mu^2*L_i^3*NumFls^2*d_i*h_i+2*mu^2*L_i^3*NumFls^2*d_i*f_i)*l)*v_i+mu^2*L_i^2*NumFls^2*brPoint^2*d_i^2*
                    l^2*slope^2-2*mu^2*L_i^2*NumFls^2*brPoint*d_i^2*l^2*slope+mu^2*L_i^2*NumFls^2*d_i^2*l^2)*beta_i^2+
                    (((-2*mu^2*L_i^3*L_max_i*NumFls^2*brPoint^2*l^2*p_i^2-4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*brPoint^2*l^2*p_i-2*
                    mu^2*H_i^2*L_i*L_max_i*NumFls^2*brPoint^2*l^2)*slope^2+(4*mu^2*L_i^3*L_max_i*NumFls^2*brPoint*l^2*p_i^2+(-8*mu*
                    L_i^3*L_max_i*NumFls*l^2*n_i+8*mu^2*H_i*L_i^2*L_max_i*NumFls^2*brPoint*l^2+(4*mu^2*L_i^3*L_max_i*NumFls^2*
                    brPoint*h_i+4*mu^2*L_i^3*L_max_i*NumFls^2*brPoint*f_i)*l)*p_i-8*mu*H_i*L_i^2*L_max_i*NumFls*l^2*n_i+4*mu^2*
                    H_i^2*L_i*L_max_i*NumFls^2*brPoint*l^2+(4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*brPoint*h_i+4*mu^2*H_i*L_i^2*
                    L_max_i*NumFls^2*brPoint*f_i)*l)*slope-2*mu^2*L_i^3*L_max_i*NumFls^2*l^2*p_i^2+((-4*mu^2*L_i^3*L_max_i*
                    NumFls^2*h_i-4*mu^2*L_i^3*L_max_i*NumFls^2*f_i)*l-4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*l^2)*p_i-2*mu^2*H_i^2*
                    L_i*L_max_i*NumFls^2*l^2+(-4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*h_i-4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*f_i)*l-2*
                    mu^2*L_i^3*L_max_i*NumFls^2*h_i^2-4*mu^2*L_i^3*L_max_i*NumFls^2*f_i*h_i-2*mu^2*L_i^3*L_max_i*NumFls^2*
                    f_i^2)*v_i^2+((-6*mu^2*L_i^2*L_max_i*NumFls^2*brPoint^2*d_i*l^2*p_i-6*mu^2*H_i*L_i*L_max_i*NumFls^2*
                    brPoint^2*d_i*l^2)*slope^2+(12*mu^2*L_i^2*L_max_i*NumFls^2*brPoint*d_i*l^2*p_i-12*mu*L_i^2*L_max_i*NumFls*
                    d_i*l^2*n_i+12*mu^2*H_i*L_i*L_max_i*NumFls^2*brPoint*d_i*l^2+(6*mu^2*L_i^2*L_max_i*NumFls^2*brPoint*d_i*
                    h_i+6*mu^2*L_i^2*L_max_i*NumFls^2*brPoint*d_i*f_i)*l)*slope-6*mu^2*L_i^2*L_max_i*NumFls^2*d_i*l^2*p_i-6*
                    mu^2*H_i*L_i*L_max_i*NumFls^2*d_i*l^2+(-6*mu^2*L_i^2*L_max_i*NumFls^2*d_i*h_i-6*mu^2*L_i^2*L_max_i*
                    NumFls^2*d_i*f_i)*l)*v_i-4*mu^2*L_i*L_max_i*NumFls^2*brPoint^2*d_i^2*l^2*slope^2+8*mu^2*L_i*L_max_i*
                    NumFls^2*brPoint*d_i^2*l^2*slope-4*mu^2*L_i*L_max_i*NumFls^2*d_i^2*l^2)*beta_i+((mu^2*L_i^2*L_max_i^2*
                    NumFls^2*brPoint^2*l^2*p_i^2+2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*brPoint^2*l^2*p_i+mu^2*H_i^2*L_max_i^2*
                    NumFls^2*brPoint^2*l^2)*slope^2+(-2*mu^2*L_i^2*L_max_i^2*NumFls^2*brPoint*l^2*p_i^2+(4*mu*L_i^2*L_max_i^2*
                    NumFls*l^2*n_i-4*mu^2*H_i*L_i*L_max_i^2*NumFls^2*brPoint*l^2+(-2*mu^2*L_i^2*L_max_i^2*NumFls^2*brPoint*h_i-2*
                    mu^2*L_i^2*L_max_i^2*NumFls^2*brPoint*f_i)*l)*p_i+4*mu*H_i*L_i*L_max_i^2*NumFls*l^2*n_i-2*mu^2*H_i^2*
                    L_max_i^2*NumFls^2*brPoint*l^2+(-2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*brPoint*h_i-2*mu^2*H_i*L_i*L_max_i^2*
                    NumFls^2*brPoint*f_i)*l)*slope+mu^2*L_i^2*L_max_i^2*NumFls^2*l^2*p_i^2+(2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*
                    l^2+(2*mu^2*L_i^2*L_max_i^2*NumFls^2*h_i+2*mu^2*L_i^2*L_max_i^2*NumFls^2*f_i)*l)*p_i+mu^2*H_i^2*
                    L_max_i^2*NumFls^2*l^2+(2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*h_i+2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*f_i)*
                    l+mu^2*L_i^2*L_max_i^2*NumFls^2*h_i^2+2*mu^2*L_i^2*L_max_i^2*NumFls^2*f_i*h_i+mu^2*L_i^2*L_max_i^2*NumFls^2*
                    f_i^2)*v_i^2+((4*mu^2*L_i*L_max_i^2*NumFls^2*brPoint^2*d_i*l^2*p_i+4*mu^2*H_i*L_max_i^2*NumFls^2*brPoint^2*
                    d_i*l^2)*slope^2+(-8*mu^2*L_i*L_max_i^2*NumFls^2*brPoint*d_i*l^2*p_i+8*mu*L_i*L_max_i^2*NumFls*d_i*l^2*
                    n_i-8*mu^2*H_i*L_max_i^2*NumFls^2*brPoint*d_i*l^2+(-4*mu^2*L_i*L_max_i^2*NumFls^2*brPoint*d_i*h_i-4*mu^2*
                    L_i*L_max_i^2*NumFls^2*brPoint*d_i*f_i)*l)*slope+4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*l^2*p_i+4*mu^2*H_i*
                    L_max_i^2*NumFls^2*d_i*l^2+(4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*h_i+4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*f_i)*
                    l)*v_i+4*mu^2*L_max_i^2*NumFls^2*brPoint^2*d_i^2*l^2*slope^2-8*mu^2*L_max_i^2*NumFls^2*brPoint*d_i^2*l^2*
                    slope+4*mu^2*L_max_i^2*NumFls^2*d_i^2*l^2)+(((mu*L_i^2*NumFls*brPoint*l*p_i+mu*H_i*L_i*NumFls*brPoint*l)*
                    slope-mu*L_i^2*NumFls*l*p_i-mu*H_i*L_i*NumFls*l+mu*L_i^2*NumFls*h_i+mu*L_i^2*NumFls*f_i)*v_i+mu*L_i*
                    NumFls*brPoint*d_i*l*slope-mu*L_i*NumFls*d_i*l)*beta_i+((-mu*L_i*L_max_i*NumFls*brPoint*l*p_i-mu*H_i*
                    L_max_i*NumFls*brPoint*l)*slope+mu*L_i*L_max_i*NumFls*l*p_i+mu*H_i*L_max_i*NumFls*l-mu*L_i*L_max_i*NumFls*
                    h_i-mu*L_i*L_max_i*NumFls*f_i)*v_i-2*mu*L_max_i*NumFls*brPoint*d_i*l*slope+2*mu*L_max_i*NumFls*d_i*l)/(((2*
                    mu*L_i^2*NumFls*l*p_i+2*mu*H_i*L_i*NumFls*l)*v_i+2*mu*L_i*NumFls*d_i*l)*beta_i+(-2*mu*L_i*L_max_i*NumFls*l*
                    p_i-2*mu*H_i*L_max_i*NumFls*l)*v_i-4*mu*L_max_i*NumFls*d_i*l)
                   Sval=ifelse(Sval>1,1,Sval)
                 }
                 Sval
                },
             stop(paste0('Foraging type: "',forageType,'" not found.'))
    )
  }

  #If S is greater than 1, sets it back to 1 (flowers can't give more than their max).
  if(S>1) S=1

  #Calculate currency
  # Multi-nest version
  # #Vector of currencies for cell, using getCurr to return appropriate nest-specific currency
  # currs=mapply(getCurr,whatCurr=whatCurr_i, L=L_i,L_max=L_max_i, e=e, d=d_i, v=v_i, h=h_i,
  #              f=f_i, l=l, p_i=p_i, c_i=c_i, c_f=c_f, H=H_i, beta=beta_i,S)
  #Single-nest version
  currs=getCurr(whatCurr_i,L_i,L_max_i,e,d_i,v_i,h_i,f_i,l,p_i,c_i,c_f,H_i,beta_i,S)

  #If sumAll is T, add all currencies together and return (useful for optimization of L).
  #Else, returns vector of currencies and S value for that cell (useful for generally returning currency)
  if(sumAll) return (sum(currs)) else return(setNames(c(currs,S),c(whatCurr_i,'S')))
}
