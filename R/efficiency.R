#Efficiency = (Gains - Foraging Loss - Travel Loss - Hive Loss)/Loading Loss + Travel Loss + Hive Loss
efficiency=function(L,L_max,e,d,v,h,l,p_i,c_i,c_f,H,S){
  Gains=L*e
  FlightLoss=(d*c_f/v)*(2+alpha(c_f,L_max,e)*(L/L_max))
  #ForagingLoss=(L*(c_i*(h+S*l*p_i)+(c_f*alpha(c_f,L_max,e)/2*L_max)*(L+S*l)))/(S*l)#Old formula
  ForagingLoss=L*(2*L_max*S*c_i*l*p_i+L*alpha(c_f,L_max,e)*c_f*l^2+S*alpha(c_f,L_max,e)*c_f*l+2*L_max*c_i*h)/(2*L_max*S*l)
  HiveLoss=c_i*H
  return((Gains-FlightLoss-ForagingLoss-HiveLoss)/(FlightLoss+ForagingLoss+HiveLoss))
}
