#Nonlinear system for random_nonsat, singleNN_nonsat, and multiNN_sat using a,
#b, c parameters. Derived from Possingham 1988 Eq 6.
nonlinearS=function(S,A,B,C,mu,l,NumFls,L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i){
  abs(((A*C*((L_i^2*l*n_i*v_i*beta_i-L_i*l*n_i*v_i)/(((mu*L_i^2*NumFls*S*l*p_i+mu*H_i*L_i*NumFls*S*l+mu*L_i^2*NumFls*h_i+mu*L_i^2*NumFls*f_i)*v_i+mu*L_i*NumFls*S*d_i*l)*beta_i+(-mu*L_i*NumFls*S*l*p_i-mu*H_i*NumFls*S*l-mu*L_i*NumFls*h_i-mu*L_i*NumFls*f_i)*v_i-2*mu*NumFls*S*d_i*l))^(B)+C+1)/(A*((L_i^2*l*n_i*v_i*beta_i-L_i*l*n_i*v_i)/(((mu*L_i^2*NumFls*S*l*p_i+mu*H_i*L_i*NumFls*S*l+mu*L_i^2*NumFls*h_i+mu*L_i^2*NumFls*f_i)*v_i+mu*L_i*NumFls*S*d_i*l)*beta_i+(-mu*L_i*NumFls*S*l*p_i-mu*H_i*NumFls*S*l-mu*L_i*NumFls*h_i-mu*L_i*NumFls*f_i)*v_i-2*mu*NumFls*S*d_i*l))^(B)+1))-S)
}
