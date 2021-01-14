#'Get currency (or currencies)
#'
#'Function to calculate currency (or currencies) in a given cell. Called by \code{curr} to return nest-specific currency.
#'
#'@param whatCurr Currency to use. Must be either \code{"rat"} (net rate) or \code{"eff"} (efficiency)
#'@param L Load size (\eqn{\muL})
#'@param L_max Maximum load size (\eqn{\muL})
#'@param e Energetic value of nectar (J/\eqn{\muL})
#'@param d Distance from hive (m)
#'@param v Unloaded flight speed from hive (m/s)
#'@param h Handling time per flower (s)
#'@param f Flight time between flowers (s)
#'@param l Maximum standing crop per flower (\eqn{\muL})
#'@param p_i Licking speed for nectar (\eqn{\muL/s})
#'@param c_i Cost of non-flying behaviour (J/s)
#'@param c_f Cost of flight (J/s)
#'@param H Time spent inside hive (s)
#'@param alphaVal Increase in metabolic rate with load (J/(s*\eqn{\muL}))
#'@param betaVal Decrease in flight speed with load (m/s\eqn{\muL})
#'@param S Competition term (0-1)
#'
#'
#'@return Currency, or vector of currencies
#'
#'@examples
#'#Parameters
#'params=list(L=59.5,
#'L_max=59.5,
#'e=14,
#'d=100,
#'v=7.8,
#'h=1.5,
#'f=0.86,
#'l=1,
#'p_i=1,
#'c_i=0.0042,
#'c_f=0.05,
#'H=100,
#'alphaVal=0.013,
#'betaVal=0.102/59.5)
#'
#'with(params,getCurr(whatCurr='rat',L,L_max,e,d,v,h,f,
#'            l,p_i,c_i,c_f,H,betaVal,S=1))
#'with(params,getCurr(whatCurr='rat',L,L_max,e,d,v,h,f,
#'            l,p_i,c_i,c_f,H,betaVal,S=0.5))
#'

getCurr=function(whatCurr,L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S){
  if(is.na(whatCurr)|is.null(whatCurr)) stop('Currency not defined')
  switch(whatCurr,
         rat=netRate(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S),
         eff=efficiency(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S),
         stop(whatCurr,' Currency not defined'))
}
