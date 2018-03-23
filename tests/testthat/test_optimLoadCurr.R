context('Currency-load optimization function (optimLoadCurr)')

nests<-list(nest1=list(xloc=1,yloc=1,n=382,whatCurr="eff",sol=FALSE,
             eps=0,L_max=59.5,v=7.8,beta=0.102,p_i=1,
             h=1.5,c_f=0.05,c_i=0.0042,H=100,d=10,
             L=59.5,curr=0))
world<-list(mu=8.33e-05,flDens=520,e=14.3,l=1,
            f=0.86,cellSize=10,forageType='random',S=1)
u<-1
scenario<-list(nests=nests,world=world)

test_that('Currency calculations work properly',{

  #Random foraging

  #Efficiency maximizers
  test1<-optimLoadCurr(u,scenario)

  expect_equal(test1$optimL,4.960216,tol=1e-04)
  expect_equal(test1$optimCurr,40.24228,tol=1e-04)
  expect_equal(test1$S,0.03390579,tol=1e-04)

  #Rate maximizers
  scenario$nests[[1]]$whatCurr='rat'
  test1b<-optimLoadCurr(u,scenario)

  #Despite having a theoretical load size of L_max, rate-maximizers have similar
  #optimal Load size to efficiency-maximizers (at least a high forager numbers)
  expect_equal(test1b$optimL,8.476394,tol=1e-04)
  expect_equal(test1b$optimCurr,0.1532347,tol=1e-04)
  expect_equal(test1b$S,0.03031699,tol=1e-04)

  #Omniscient foraging style:

  #Efficiency maximizers
  scenario$world$forageType='omniscient'
  scenario$nests[[1]]$whatCurr='eff'

  test1<-optimLoadCurr(u,scenario)

  expect_equal(test1$optimL,4.915216,tol=1e-04)
  expect_equal(test1$optimCurr,41.50074,tol=1e-04)
  expect_equal(test1$S,0.03558491,tol=1e-04)

  #Rate maximizers
  scenario$nests[[1]]$whatCurr='rat'
  test1b<-optimLoadCurr(u,scenario)

  expect_equal(test1b$optimL,4.91522,tol=1e-04)
  expect_equal(test1b$optimCurr,0.1583363,tol=1e-04)
  expect_equal(test1b$S,0.0355849,tol=1e-04)

})

test_that('Exception handling works properly',{
  #Empty nests in non-worthless patch
  scenario$nests[[1]]$n=0 #Sets forager number to 0
  test2<-optimLoadCurr(u,scenario)
  expect_equivalent(test2$optimL,NA) #L should be NA
  expect_equivalent(test2$optimCurr,0) #curr should be 0
  expect_identical(test2$S,1) #S should be 1
  scenario$nests[[1]]$n=382 #Resets forager number

  #Worthless patch
  scenario$world$mu=0 #Sets per flower production to 0
  test3<-optimLoadCurr(u,scenario)
  expect_equivalent(test3$optimL,NA) #L should be NA
  expect_equivalent(test3$optimCurr,-1) #curr should be -1 (or -Inf for rate maximizers)
  expect_identical(test3$S,NA) #S should be NA
  scenario$world$mu=8.33e-05

  #Forager arguments NA
  scenario$nests[[1]]$whatCurr=NA
  expect_error(optimLoadCurr(u,scenario))

  #Forager arguments missing
  scenario$nests[[1]]$whatCurr=NULL
  expect_error(optimLoadCurr(u,scenario))
  scenario$nests[[1]]$whatCurr='eff'

  #Too many cells provided
  expect_error(optimLoadCurr(u=c(1,2),scenario))
})

