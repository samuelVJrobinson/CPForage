#'Change in loaded flight cost (\eqn{\alpha})
#'
#'Returns \eqn{\alpha}, which is the proportion increase in completely loaded cost of flight
#'
#'@param c_f Cost of (unloaded) flight (J/s)
#'@param L_max Maximum load size (\eqn{\muL})
#'@param e_i Energetic value of load (\eqn{J/\muL})
#'@param alphaVal Proportion cost increase (default: 5e-05 from Schmid-Hempel et al 1985)
#'@return \eqn{\alpha} value: units \eqn{J/\muL}
#'@examples
#'alpha(0.05,59.5,14.35)
#'
alpha=function(c_f,L_max,e_i,alphaVal=5e-05) (c_f+L_max*e_i*alphaVal)/c_f-1

#NOTE: I think this should be (c_f+L_max*e_i*alphaVal)/c_f. Alternatively, this entire expression could be moved outside of the function
