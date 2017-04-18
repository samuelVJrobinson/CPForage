context('Scaling function')

test_that("scalFun is between 0 and 1",{
  mu=0.3/3600
  l=1
  NumFls=520*(10^2)
  L_i=59.5
  L_max_i=59.5
  h_i=1.5
  f_i=0.86
  p_i=1
  d_i=100
  v_i=7.8
  beta_i=0.102
  H_i=100

  expect_equal(round(scalFun(mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100, #Patch-level
          h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=T),5),round(0.1183501,5))
  expect_equal(round(scalFun(mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100, #Flower-level
                       h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=F),5),round(0.1042959,5))
  expect_equal(round(scalFun(mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100, #Patch-level with 2 nests
                       h_i=h_i,p_i=p_i,f_i=f_i,d_i=c(d_i,200),v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=T),5),round(0.05517226,5))
  expect_equal(round(scalFun(mu=mu,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100, #Flower-level with 2 nests
                       h_i=h_i,p_i=p_i,f_i=f_i,d_i=c(d_i,200),v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=F),5),round(0.2417517,5))

  expect_error(scalFun(mu=NA,l=l,NumFls=NumFls,L_i=L_i,L_max_i=L_max_i,n_i=100, #Flower-level with 2 nests
                       h_i=h_i,p_i=p_i,f_i=f_i,d_i=c(d_i,200),v_i=v_i,beta_i=beta_i,H_i=H_i,patchLev=F))

})
