context('Summed currency selection (curr)')

params=list(L_i=59.5,
            L_max_i=59.5,
            n_i=10,
            h_i=1.5,
            p_i=1,
            f_i=0.86,
            d_i=100,
            v_i=7.8,
            beta_i=0.102,
            H_i=100,
            c_f=0.05,
            c_i=0.0042,
            mu=0.3/3600,
            l=1,
            e=14.35,
            NumFls=520*(10^2))

test_that('Summed currency works properly',{
  #Efficiency, using flower-level competition
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,patchLev=F)),166.8337,tol=1e-4)
  #Netrate, using flower-level competition
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                      c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=T,patchLev=F)),2.0684,tol=1e-4)
  #Efficiency, using patch-level competition
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                      c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,patchLev=T)),201.1507,tol=1e-4)
  #Netrate, using patch-level competition
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                      c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=T,patchLev=T)),2.5347,tol=1e-4)
  #Netrate, using patch-level competition, with 2 nests
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,c(d_i,300),v_i,beta_i,H_i,
                         c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=T,patchLev=T)),4.6995,tol=1e-4)

})

test_that('Individual currency and S-values work properly',{
  #Efficiency with 1 nest
  expect_equal(length(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                       c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,patchLev=F))),2)
  #Efficiency and rate with 2 nests
  expect_equal(length(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,c(d_i,200),v_i,beta_i,H_i,
                                       c_i,c_f,mu,l,e,NumFls,whatCurr_i=c('eff','rat'),sumAll=F,patchLev=F))),3)
})



