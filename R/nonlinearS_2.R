#Nonlinear system for multiNN_sat competion. Derived from Possingham 1988 Eq 6.
nonlinearS_2=function(S,a0,a1,b1,mu,l,NumFls,L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i){
  abs(-((L_i*a1*l)/(((L_i^2*l*n_i*v_i*beta_i-L_i*l*n_i*v_i)/(((mu*L_i^2*NumFls*S*l*p_i+mu*H_i*L_i*NumFls*S*l+mu*L_i^2*NumFls*h_i+mu*L_i^2*NumFls*f_i)*v_i+mu*L_i*NumFls*S*d_i*l)*beta_i+(-mu*L_i*NumFls*S*l*p_i-mu*H_i*NumFls*S*l-mu*L_i*NumFls*h_i-mu*L_i*NumFls*f_i)*v_i-2*mu*NumFls*S*d_i*l))^((L_i*b1*l-S)/S))-1)/((a0)/(((L_i^2*l*n_i*v_i*beta_i-L_i*l*n_i*v_i)/(((mu*L_i^2*NumFls*S*l*p_i+mu*H_i*L_i*NumFls*S*l+mu*L_i^2*NumFls*h_i+mu*L_i^2*NumFls*f_i)*v_i+mu*L_i*NumFls*S*d_i*l)*beta_i+(-mu*L_i*NumFls*S*l*p_i-mu*H_i*NumFls*S*l-mu*L_i*NumFls*h_i-mu*L_i*NumFls*f_i)*v_i-2*mu*NumFls*S*d_i*l))^((L_i*b1*l-S)/S))+1)-S)
}
