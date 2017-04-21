#'Change in loaded flight cost (\eqn{\alpha})
#'
#'Returns \eqn{\alpha}, which is the proportion increase in completely loaded cost of flight (5e-05 value from Schmid-Hempel et al 1985)
#'
#'@param c_f Cost of (unloaded) flight (J/s)
#'@param L_max Maximum load size (\eqn{\muL})
#'@param e_i Energetic value of load (\eqn{J/\muL})
#'@return \eqn{\alpha} value
#'@examples
#'alpha(0.05,59.5,14.35)
#'
alpha=function(c_f,L_max,e_i) (c_f+L_max*e_i*5e-05)/c_f-1
