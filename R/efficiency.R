#'Efficiency
#'
#'Returns Efficiency (ratio of profits/losses) for a forager in a given cell.
#'
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
#'@param S Competition term (0-1)
#' @param alphaVal Increase in metabolic rate with load (J/(s*\eqn{\muL}))
#'@param betaVal Reduction in flight speed with increase in load
#'
#'@return Efficiency (dimensionless). \eqn{Efficiency = \frac{Gains - Foraging
#'  Loss - Travel Loss - Hive Loss}{Loading Loss + Travel Loss + Hive Loss}}
#'  Called by \code{curr_i}.
#'@examples
#'efficiency(L=50,L_max=59.5,e=14.35,d=100,v=7.8,
#'  h=1.5,f=0.86,l=1,p_i=1,c_i=0.0042,c_f=0.05,H=100,S=0.5,
#'  alphaVal=0.013,betaVal=0.102/59.5)
#'

efficiency <- function(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S){

  Gains <- L*e #Gains within patch

  OutboundLoss <- c_f*d/v #Cost of traveling to patch from hive
  InboundLoss <- ((L*alphaVal+1)*(L*betaVal+1)*c_f*d)/v #Cost of traveling back to hive
  FlightLoss <- OutboundLoss+InboundLoss #Total flight costs

  ForageLossHandling <- (L*c_i*(S*l*p_i+h))/(S*l) #Energy required to extract necter

  nFlwVis <- L/(S*l) #Number of flowers visited in a trip

  if(nFlwVis<1){
    ForageLossFlying <- 0  #Only 1 flower visited, so no intra-patch movement needed
  } else {
    #Energetic losses while flying from flower-to-flower
    ForageLossFlying <- c_f*f*(
      (2*S^2*alphaVal*betaVal*l^2*(nFlwVis-1)^3+3*S^2*alphaVal*betaVal*l^2*
         (nFlwVis-1)^2+S^2*alphaVal*betaVal*l^2*(nFlwVis-1))/6+
           (S*betaVal*l*(nFlwVis-1)^2+S*betaVal*l*(nFlwVis-1))/2+
           (S*alphaVal*l*(nFlwVis-1)^2+S*alphaVal*l*(nFlwVis-1))/2+
           nFlwVis-1)
  }

  ForagingLoss <- ForageLossHandling+ForageLossFlying #Total foraging loss
  HiveLoss <- c_i*H #Loss within hive

  #BS checking
  if(any(c(Gains,FlightLoss,ForagingLoss,HiveLoss)<0)){ #Gains/losses can't be negative
    stop(paste0('Efficiency calcs failed:\nGains=',Gains,'\nFlightLoss=',FlightLoss,
                '\nForagingLoss=',ForagingLoss,'\nHiveLoss=',HiveLoss))
  }

  Efficiency <- (Gains-FlightLoss-ForagingLoss-HiveLoss)/(FlightLoss+ForagingLoss+HiveLoss)

  return(Efficiency)
}



