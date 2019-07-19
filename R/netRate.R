#' Net rate
#'
#' Helper function for getCurr. Returns Net rate of Energy Gain (J/s) for a
#' forager in a given cell.
#'
#' @param L Load size (\eqn{\muL})
#' @param L_max Maximum load size (\eqn{\muL})
#' @param e Energetic value of nectar (J/\eqn{\muL})
#' @param d Distance from hive (m)
#' @param v Unloaded flight speed from hive (m/s)
#' @param h Handling time per flower (s)
#' @param f Travel time between flowers (s)
#' @param l Maximum standing crop per flower (\eqn{\muL})
#' @param p_i Licking speed for nectar (\eqn{\muL/s})
#' @param c_i Cost of non-flying behaviour (J/s)
#' @param c_f Cost of flight (J/s)
#' @param H Time spent inside hive (s)
#' @param S Competition term (0-1)
#' @param beta Reduction in flight speed with increase in load
#' @param alpha Increase in flight cost with increase in load (5e-5 by default)
#'
#' @return Net rate (J/s). \eqn{Net rate = \frac{Gains - Foraging
#'  Loss - Travel Loss - Hive Loss}{Loading Time + Travel Time + Hive Time}}.
#'  Called by \code{curr_i}.
#' @export
#'
#' @examples
#'netRate(L=50,L_max=50.5,e=14.35,d=100,v=7.8,
#'  h=1.5,f=0.86,l=1,p_i=1,c_i=0.0042,c_f=0.05,H=100,beta=0.102,S=0.5)
netRate=function(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S,alpha=5e-05){
  #Rate=(Gains-PatchLoss-Travel Loss - Hive Loss)/(Travel Time + Foraging Time + Hive Time)

  #Rate maximizers "should" always take largest load, unless S is high (at that point netrate(small load)=netrate(large load))
  Gains <- L*e

  OutboundFlightLoss <- c_f*d/v #Traveling from hive to patch
  InboundFlightLoss <- (c_f+L*e*alpha)*d/v #Traveling from patch back to hive
  FlightLoss <- OutboundFlightLoss+InboundFlightLoss #Total flight costs
  FlightTime <- (d*(L*beta-2*L_max))/(v*(L*beta-L_max)) #Time taken to fly to and from patch

  PatchLossHandling <- L*c_i*(S*l*p_i+h)/S*l #Total cost of handling flowers in a patch
  PatchTimeHandling <- L*(S*l*p_i+h)/S*l #Time taken to handle flowers in a patch

  if(L/S*l<2){ #If less than 2 flowers are visited
    PatchLossFlying <- 0  #Only 1 flower visited, so no intra-patch movement needed
    PatchTimeFlying <- 0
  } else {
    PatchLossFlying <- f*((S*e*(L/(S*l)-1)*l*alpha+S*e*(L/(S*l)-1)^2*l*alpha)/2+c_f*(L/(S*l)-1)) #Energy loss while flying between flowers

    #Linear approximation of summed flight time. Not exact, but faster than summation.
    y_upr <- f*(1/(1-S*(L/S*l-1)*l*beta/L_max))
    y_lwr <- f*(1/(1-S*l*beta/L_max))
    PatchTimeFlying <- y_lwr*(L/S*l-2) + (y_upr-y_lwr)*(L/S*l-2)*0.5
  }

  PatchLoss <- PatchLossHandling+PatchLossFlying #Total loss during foraging
  PatchTime <- PatchTimeHandling+PatchTimeFlying #Total time taken to forage
  HiveLoss <- c_i*H #Energy loss while unloading in hive
  HiveTime <- H #Time taken to unload in hive
  return((Gains-FlightLoss-PatchLoss-HiveLoss)/(FlightTime+PatchTime+HiveTime))
}
