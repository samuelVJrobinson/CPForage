#' Efficiency (ratio of profits/losses)
#'
#' Calculates efficiency for a forager in a given cell.
#' \deqn{Efficiency = \frac{Gains - Foraging Loss - Travel Loss - Hive Loss}{Loading Loss + Travel Loss + Hive Loss}}
#' Called by \code{curr_i}.
#'
#'@param L Load size (\eqn{\muL})
#'@param L_max Maximum load size (\eqn{\muL})
#'@param e Energetic value of nectar (J/\eqn{\muL})
#'@param d Distance from hive (m)
#'@param v Unloaded flight speed from hive (m/s)
#'@param h Handling time per flower (s)
#'@param l Maximum standing crop per flower (\eqn{\muL})
#'@param p_i Licking speed for nectar (\eqn{\muL/s})
#'@param c_i Cost of non-flying behaviour (J/s)
#'@param c_f Cost of flight (J/s)
#'@param H Time spent inside hive (s)
#'@param S Competition term (0-1)
#'
#'@return Efficiency (dimensionless)
#'@examples
#'efficiency(L=50,L_max=50.5,e=14.35,d=100,v=7.8,
#'  h=1.5,l=1,p_i=1,c_i=0.0042,c_f=0.05,H=100,S=0.5)

efficiency=function(L,L_max,e,d,v,h,l,p_i,c_i,c_f,H,S){
  Gains=L*e
  FlightLoss=(d*c_f/v)*(2+alpha(c_f,L_max,e)*(L/L_max))
  #ForagingLoss=(L*(c_i*(h+S*l*p_i)+(c_f*alpha(c_f,L_max,e)/2*L_max)*(L+S*l)))/(S*l)#Old formula
  ForagingLoss=L*(2*L_max*S*c_i*l*p_i+L*alpha(c_f,L_max,e)*c_f*l^2+S*alpha(c_f,L_max,e)*c_f*l+2*L_max*c_i*h)/(2*L_max*S*l)
  HiveLoss=c_i*H
  return((Gains-FlightLoss-ForagingLoss-HiveLoss)/(FlightLoss+ForagingLoss+HiveLoss))
}
