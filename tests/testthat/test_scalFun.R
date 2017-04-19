context('Scaling function (scalFun)')

test_that("scalFun is between 0 and 1",{
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

  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100, #Patch-level
          h_i,p_i,f_i,d_i,v_i,beta_i,H_i,patchLev=T)),0.1183501,tol=1e-4)
  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100, #Flower-level
                       h_i,p_i,f_i,d_i,v_i,beta_i,H_i,patchLev=F)),0.1042959,tol=1e-4)
  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100, #Patch-level with 2 nests
                                   h_i,p_i,f_i,d_i=c(d_i,200),v_i,beta_i,H_i,patchLev=T)),0.05517226,tol=1e-4)
  expect_equal(with(params,scalFun(mu,l,NumFls,L_i,L_max_i,n_i=100, #Flower-level with 2 nests
                                   h_i,p_i,f_i,d_i=c(d_i,200),v_i,beta_i,H_i,patchLev=F)),0.2417517,tol=1e-4)
  })

test_that('scalFun error handling works',{
  #NA value for nectar production
  expect_error(with(params,scalFun(mu=NA,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100,
                       h_i=h_i,p_i=p_i,f_i=f_i,d_i=c(d_i,200),v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=F)))

  #scalFun "should" work for slightly negative L_i values (used by optimizer)
  expect_equal(with(params,scalFun(mu=mu,l=l,NumFls=NumFls,L_i=-0.1,L_max_i=L_max_i,n_i=100,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=F)),1)

  #However, it should not work for any other values being negative
  expect_error(with(params,scalFun(mu=mu,l=l,NumFls=-100,L_i=L_i,L_max_i=L_max_i,n_i=100,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=F)))

  expect_error(with(params,scalFun(mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=-100,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=F)))

})
