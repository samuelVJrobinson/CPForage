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
#' @param l Maximum standing crop per flower (\eqn{\muL})
#' @param p_i Licking speed for nectar (\eqn{\muL/s})
#' @param c_i Cost of non-flying behaviour (J/s)
#' @param c_f Cost of flight (J/s)
#' @param H Time spent inside hive (s)
#' @param S Competition term (0-1)
#' @param beta Reduction in flight speed with increase in load
#'
#' @return Net rate (J/s). \eqn{Net rate = \frac{Gains - Foraging
#'  Loss - Travel Loss - Hive Loss}{Loading Time + Travel Time + Hive Time}}.
#'  Called by \code{curr_i}.
#' @export
#'
#' @examples
#'netRate(L=50,L_max=50.5,e=14.35,d=100,v=7.8,
#'  h=1.5,l=1,p_i=1,c_i=0.0042,c_f=0.05,H=100,beta=0.102,S=0.5)
netRate=function(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S){
  #Rate=(Gains-ForagingLoss-Travel Loss - Hive Loss)/(Travel Time + Foraging Time + Hive Time)

  #Rate maximizers "should" always take largest load, unless S is high (at that point netrate(small load)=netrate(large load))
  Gains=L*e
  FlightLoss= (c_f*d*(L*alpha(c_f,L_max,e)+L_max))/(L_max*v)
  FlightTime=(d*(L*beta-2*L_max))/(v*(L*beta-L_max))

  ForageLossHandling=(L*c_i*(S*l*p_i+h))/(S*l)
  ForageTimeHandling=L*(S*l*p_i+h)/S*l

  if(L/S*l<2){ #If less than 2 flowers are visited
    ForageLossFlying = ForageTimeFlying = 0  #Only 1 flower visited, so no intra-patch movement needed
  } else {
    ForageLossFlying=(S*c_f*f*(L/(S*l)+(L/(S*l)-1)^2-1)*l*alpha(c_f,L_max,e))/(2*L_max)
    ForageTimeFlying=sum(f*(1/(1-S*(1:(L/S*l-1))*l*beta/L_max))) #No simple solution; sum across terms
  }

  ForagingLoss=ForageLossHandling+ForageLossFlying
  ForagingTime=ForageTimeHandling+ForageTimeFlying
  HiveLoss=c_i*H
  HiveTime=H
  return((Gains-FlightLoss-ForagingLoss-HiveLoss)/(FlightTime+ForagingTime+HiveTime))
}
