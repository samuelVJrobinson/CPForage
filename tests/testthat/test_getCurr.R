context('Currency selection function (getCurr)')

params <- list(L=59.5,
            L_max=59.5,
            e=1.2193,
            d=100,
            v=7.8,
            h=1.5,
            f=0.86,
            l=1,
            p_i=1,
            c_i=0.0042,
            c_f=0.05,
            H=100,
            alphaVal=0.013,
            betaVal=0.102/59.5
            )

test_that('Currency calculations work properly',{
  #Net rate
  expect_equal(with(params,getCurr(whatCurr='rat',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S=1)),
               0.201,tol=1e-03) #Untouched patch
  expect_equal(with(params,getCurr(whatCurr='rat',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S=0.5)),
               0.1310567,tol=1e-03) #Patch at half value

  #Efficiency
  expect_equal(with(params,getCurr(whatCurr='eff',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S=1)),
               9.957552,tol=1e-03)
  expect_equal(with(params,getCurr(whatCurr='eff',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S=0.5)),
               5.754044,tol=1e-03)

  #Error handling - no currency defined
  expect_error(with(params,getCurr(whatCurr=NA,L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S=1)))
  expect_error(with(params,getCurr(whatCurr='blah',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,alphaVal,betaVal,S=1)))
})

