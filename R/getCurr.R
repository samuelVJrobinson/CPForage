#Function to switch between currency calculations depending on what is supplied to WHATCURR. Called by curr
# getCurr=function(whatCurr,L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S){
getCurr=function(whatCurr,L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S){
  switch(whatCurr,
         rat=netRate(L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S),
         eff=efficiency(L,L_max,e,d,v,h,l,p_i,c_i,c_f,H,S))
}
