context('Currency selection function (getCurr)')

params=list(L=59.5,
            L_max=59.5,
            e=14,
            d=100,
            v=7.8,
            h=1.5,
            f=0.86,
            l=1,
            p_i=1,
            c_i=0.0042,
            c_f=0.05,
            H=100,
            beta=0.102)

test_that('Currency calculations work properly',{
  #Net rate
  expect_equal(with(params,getCurr(whatCurr='rat',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S=1)),2.534665,tol=1e-04)
  expect_equal(with(params,getCurr(whatCurr='rat',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S=0.5)),1.769789,tol=1e-04)

  #Efficiency
  expect_equal(with(params,getCurr(whatCurr='eff',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S=1)),201.1507,tol=1e-04)
  expect_equal(with(params,getCurr(whatCurr='eff',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S=0.5)),144.258,tol=1e-04)

  #Error handling - no currency defined
  expect_error(with(params,getCurr(whatCurr=NA,L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S=1)))
  expect_error(with(params,getCurr(whatCurr='blah',L,L_max,e,d,v,h,f,l,p_i,c_i,c_f,H,beta,S=1)))
})

