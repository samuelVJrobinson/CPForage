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

  #L=L_max #Rate maximizers "should" always take largest load, unless S is high (at that point netrate(small load)=netrate(large load))
  Gains=L*e
  FlightLoss=(d*c_f/v)*(2+alpha(c_f,L_max,e)*(L/L_max))
  ForagingLoss=((L*S*c_f*f*l+L^2*c_f*f)*alpha(c_f,L_max,e)+2*L*L_max*S*c_i*l*p_i+2*L*L_max*c_i*h)/(2*L_max*S*l)
  HiveLoss=c_i*H

  FlightTime=(d/v)*((2-beta*L/L_max)/(1-beta*L/L_max))
  ForagingTime=(L*(S*l*p_i+h+f))/(S*l)
  HiveTime=H
  return((Gains-FlightLoss-ForagingLoss-HiveLoss)/(FlightTime+ForagingTime+HiveTime))
}
