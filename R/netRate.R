#Rate=(Gains-ForagingLoss-Travel Loss - Hive Loss)/(Travel Time + Foraging Time + Hive Time)
netRate=function(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S){
  #L=L_max #Rate maximizers "should" always take largest load, unless S is high (at that point netrate(small load)=netrate(large load))
  Gains=L*e
  FlightLoss=(d*c_f/v)*(2+alpha(c_f,L_max,e)*(L/L_max))
  ForagingLoss=L*(2*L_max*S*c_i*l*p_i+L*alpha(c_f,L_max,e)*c_f*l^2+S*alpha(c_f,L_max,e)*c_f*l+2*L_max*c_i*h)/(2*L_max*S*l)
  HiveLoss=c_i*H

  FlightTime=(d/v)*((2-beta*L/L_max)/(1-beta*L/L_max))
  ForagingTime=L*(h+S*l*p_i+f)/S*l
  HiveTime=H
  return((Gains-FlightLoss-ForagingLoss-HiveLoss)/(FlightTime+ForagingTime+HiveTime))
}
