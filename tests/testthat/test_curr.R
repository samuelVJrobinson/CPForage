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
  #Efficiency, using omniscient foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
              c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,forageType='omniscient')),203.9737,tol=1e-4)
  #Efficiency, using random foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
              c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,forageType='random')),163.5881,tol=1e-4)
  #Net rate, using omniscient foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
              c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=T,forageType='omniscient')),2.59821,tol=1e-4)
  #Net rate, using random foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
              c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=T,forageType='random')),2.048508,tol=1e-4)
})

test_that('Individual currency and S-values work properly',{
  #Efficiency with 1 nest
  effRes1 <- with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                   c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,forageType='random'))
  expect_equal(length(effRes1),2)
  expect_equal(effRes1[['S']],0.6167597,tol=1e-4)
  expect_equal(effRes1[['eff']],163.5881,tol=1e-4)

  #Same thing, but 200 foragers
  effRes2 <- with(params,curr(L_i,L_max_i,200,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                              c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,forageType='random'))
  expect_equal(length(effRes2),2)
  expect_equal(effRes2[['S']],0.04695388,tol=1e-4)
  expect_equal(effRes2[['eff']],21.73642420,tol=1e-4)


  #Net rate with 1 nests
  ratRes1 <- with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                              c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='random'))
  expect_equal(length(ratRes1),2)
  expect_equal(ratRes1[['S']],0.6167597,tol=1e-4)
  expect_equal(ratRes1[['rat']],2.048508,tol=1e-4)

  #Same thing, but 200 foragers
  ratRes2 <- with(params,curr(L_i,L_max_i,200,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                              c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='random'))
  expect_equal(length(ratRes2),2)
  expect_equal(ratRes2[['S']],0.04695388,tol=1e-4)
  expect_equal(ratRes2[['rat']],0.25691618,tol=1e-4)

  #NN foraging:
  #Net rate with 1 nests
  ratRes1 <- with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                              c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='random'))
  expect_equal(length(ratRes1),2)
  expect_equal(ratRes1[['S']],0.6167597,tol=1e-4)
  expect_equal(ratRes1[['rat']],2.048508,tol=1e-4)

  #Same thing, but 200 foragers
  ratRes2 <- with(params,curr(L_i,L_max_i,200,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                              c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='random'))
  expect_equal(length(ratRes2),2)
  expect_equal(ratRes2[['S']],0.04695388,tol=1e-4)
  expect_equal(ratRes2[['rat']],0.25691618,tol=1e-4)
})


test_that('Exception handling (limited - mostly in optimLoadCurr)',{
  #Load of 0
  effRes1 <- with(params,curr(0,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                              c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,forageType='random'))
  expect_equal(effRes1[['S']],1)
  expect_equal(effRes1[['eff']],-1)
})
