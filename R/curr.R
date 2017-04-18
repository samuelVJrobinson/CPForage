#Function to calculate currency(-ies) in a given cell. Called by optimLoadCurr
curr=function(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,c_i,c_f,whatCurr_i,mu,l,e,NumFls,sumAll=T,Spatch=F){
  Stemp=scalFun(mu=mu, l=l, L_i=L_i, L_max=L_max_i, n_i=n_i, #Scaling function for cell
                h_i=h_i, p_i=p_i, f_i=f_i, d_i=d_i, v_i=v_i, beta_i=beta_i, H_i=H_i,NumFls=NumFls,patchLev=Spatch)
  #Vector of currencies for cell, using getCurr to return appropriate nest-specific currency
  currs=mapply(getCurr,whatCurr=whatCurr_i, L=L_i,L_max=L_max_i, e=e, d=d_i, v=v_i, h=h_i,
               f=f_i, l=l, p_i=p_i, c_i=c_i, c_f=c_f, H=H_i, beta=beta_i,Stemp)
  #If sumAll is T, add all currencies together and return (useful for optimization of L).
  #Else, returns vector of currencies and S value for that cell (useful for generally returning currency)
  if(sumAll) return (sum(currs)) else return(setNames(c(currs,Stemp),c(names(currs),'S')))
}
