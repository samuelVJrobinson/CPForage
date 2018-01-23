context('Scaling function (scalFun)')

params=list(mu=0.3/3600,
            l=1,
            NumFls=520*(10^2),
            L_i=59.5,
            L_max_i=59.5,
            h_i=1.5,
            f_i=0.86,
            p_i=1,
            d_i=100,
            v_i=7.8,
            beta_i=0.102,
            H_i=100)

test_that("scalFun is between 0 and 1",{
  #Omniscient foraging, with non-saturating nectar production
  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100,
          h_i,p_i,f_i,d_i,v_i,beta_i,H_i,forageType='omniscient')),0.1183501,tol=1e-4)

  #Random foraging with saturating nectar production
  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100,
                       h_i,p_i,f_i,d_i,v_i,beta_i,H_i,forageType='random')),0.09459352,tol=1e-4)

  #Random foraging with non-saturating nectar production
  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100,
                      h_i,p_i,f_i,d_i,v_i,beta_i,H_i,forageType='random_nonsat')),0.1113888,tol=1e-4)

  #Single nearest-neighbour foraging with non-saturating nectar production
  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100,
                                   h_i,p_i,f_i,d_i,v_i,beta_i,H_i,forageType='singleNN_nonsat')),0.2531934,tol=1e-4)

  })

test_that('scalFun error handling works',{
  #Multiple nests aren't currently allowed
  expect_error(with(params,scalFun(mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100,
                           h_i=h_i,p_i=p_i,f_i=f_i,d_i=c(d_i,150),v_i=v_i,beta_i=beta_i,H_i=H_i,forageType='omniscient')))

  #NA value for nectar production - should crash
  expect_error(with(params,scalFun(mu=NA,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100,
                       h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,forageType='omniscient')))

  #scalFun "should" work for slightly negative L_i values (used by optimizer)
  expect_equal(with(params,scalFun(mu=mu,l=l,NumFls=NumFls,L_i=-0.1,L_max_i=L_max_i,n_i=100,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,forageType='omniscient')),1)

  #However, it should not work for any other values being negative
  expect_error(with(params,scalFun(mu=mu,l=l,NumFls=-100,L_i=L_i,L_max_i=L_max_i,n_i=100,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,forageType='omniscient')))

  expect_error(with(params,scalFun(mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=-100,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,forageType='omniscient')))

})
