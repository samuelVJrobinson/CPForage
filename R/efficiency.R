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
#'@param f Travel time between flowers (s)
#'@param l Maximum standing crop per flower (\eqn{\muL})
#'@param p_i Licking speed for nectar (\eqn{\muL/s})
#'@param c_i Cost of non-flying behaviour (J/s)
#'@param c_f Cost of flight (J/s)
#'@param H Time spent inside hive (s)
#'@param S Competition term (0-1)
#'@param alpha Increase in flight cost with increase in load (5e-5 by default)
#'
#'@return Efficiency (dimensionless). \eqn{Efficiency = \frac{Gains - Foraging
#'  Loss - Travel Loss - Hive Loss}{Loading Loss + Travel Loss + Hive Loss}}
#'  Called by \code{curr_i}.
#'@examples
#'efficiency(L=50,L_max=50.5,e=14.35,d=100,v=7.8,
#'  h=1.5,f=0.86,l=1,p_i=1,c_i=0.0042,c_f=0.05,H=100,S=0.5)

efficiency=function(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,S,alpha=5e-05){
  Gains=L*e #Energy gains in patch

  OutboundFlightLoss <- c_f*d/v #Traveling from hive to patch
  InboundFlightLoss <- (c_f+L*e*alpha)*d/v #Traveling from patch back to hive
  FlightLoss <- OutboundFlightLoss+InboundFlightLoss #Total flight costs

  PatchLossHandling <- L*c_i*(S*l*p_i+h)/S*l #Total cost of handling flowers in a patch

  if(L/S*l<2){ #If less than 2 flowers are visited
    PatchLossFlying <- 0  #Only 1 flower visited, so no intra-patch movement needed
  } else {
    PatchLossFlying <- f*((S*e*(L/(S*l)-1)*l*alpha+S*e*(L/(S*l)-1)^2*l*alpha)/2+c_f*(L/(S*l)-1)) #Energy loss while flying between flowers
  }

  PatchLoss <- PatchLossHandling+PatchLossFlying #Total energy loss while foraging
  HiveLoss <- c_i*H #Energy loss while unloading in hive
  return((Gains-FlightLoss-PatchLoss-HiveLoss)/(FlightLoss+PatchLoss+HiveLoss))
}
