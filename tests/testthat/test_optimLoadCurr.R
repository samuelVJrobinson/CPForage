context('Currency-load optimization function (optimLoadCurr)')

nests<-list(nest1=list(xloc=1,yloc=1,n=382,whatCurr="eff",sol=FALSE,
             eps=0,L_max=59.5,v=7.8,beta=0.102,p_i=1,
             h=1.5,c_f=0.05,c_i=0.0042,H=100,d=10,
             L=59.5,curr=0))
world<-list(mu=8.33e-05,flDens=520,e=14.3,l=1,
            f=0.86,cellSize=10,patchLev=FALSE,S=1)
u<-1
scenario<-list(nests=nests,world=world)

test_that('Currency calculations work properly',{
  test1<-optimLoadCurr(u,scenario)

  expect_equal(unname(test1$optimL),4.605306,tol=1e-04)
  expect_equal(names(test1$optimL),'nest1')

  expect_equal(unname(test1$optimCurr),39.4235,tol=1e-04)
  expect_equal(names(test1$optimCurr),'nest1')

  expect_equal(unname(test1$S),0.03466368,tol=1e-04)
  expect_equal(names(test1$S),'S')

  scenario$nests[[1]]$whatCurr='rat' #Rate maximizers
  test1b<-optimLoadCurr(u,scenario)

  expect_equal(unname(test1b$optimL),7.915367,tol=1e-04)

  expect_equal(unname(test1b$optimCurr),0.1530864,tol=1e-04)

  expect_equal(unname(test1b$S),0.03064256,tol=1e-04)
  expect_equal(names(test1b$S),'S')

  #Rate vs Rate
  scenario$nests$nest2=scenario$nests$nest1 #Add 2nd nest
  test1c<-optimLoadCurr(u,scenario)
  expect_equal(unname(test1c$optimL),c(16.04992,16.04992 ),tol=1e-04)
  expect_equal(names(test1c$optimL),c('nest1','nest2'))
  expect_equal(unname(test1c$optimCurr),c(0.3092051,0.3092051),tol=1e-04)
  expect_equal(names(test1c$optimCurr),c('nest1','nest2'))
  expect_equal(unname(test1c$S),0.06194489,tol=1e-04)
  expect_equal(names(test1c$S),'S')

  #Rate vs Eff
  scenario$nests$nest2$whatCurr='eff' #Changes 2nd nest to Eff
  test1d<-optimLoadCurr(u,scenario)
  expect_equal(unname(test1d$optimL),c(0.07769018,38.93359533),tol=1e-04)
  expect_equal(unname(test1d$optimCurr),c(5.464097e-03,3.585243e+02),tol=1e-04)
  expect_equal(unname(test1d$S),0.9999339,tol=1e-04)
  #NOTE: THE RESULTS FROM THIS ARE WEIRD, BUT SOCIAL EFFIENCY MAXIMIZERS ALWAYS WIN AGAINST SOCIAL RATE MAXIMIZERS

  #Rate vs Eff (solitary) - acts the same as social because of "social" optimization criteria (maximizing output of cell)
  scenario$nests$nest1$sol=T #Changes both nests to solitaries
  scenario$nests$nest2$sol=T
  test1e<-optimLoadCurr(u,scenario)
  expect_equal(unname(test1e$optimL),c(0.07769018,38.93359533),tol=1e-04)
  expect_equal(unname(test1e$optimCurr),c(5.464097e-03,3.585243e+02),tol=1e-04)
  expect_equal(unname(test1e$S),0.9999339,tol=1e-04)
  scenario$nests$nest2=NULL
})

test_that('Exception handling works properly',{
  #Empty nests in non-worthless patch
  scenario$nests[[1]]$n=0
  test2<-optimLoadCurr(u,scenario)
  expect_identical(unname(test2$optimL),0) #L should be 0
  expect_identical(unname(test2$optimCurr),0) #curr should be 0
  expect_identical(unname(test2$S),1) #S should be 1
  scenario$nests[[1]]$n=382

  #Worthless patch
  scenario$world$mu=0
  test3<-optimLoadCurr(u,scenario)
  expect_identical(unname(test3$optimL),NA) #L should be NA
  expect_identical(unname(test3$optimCurr),-Inf) #curr should be -Inf (or -1 for eff maximizers)
  expect_identical(unname(test3$S),NA) #S should be NA
  world$mu=8.33e-05

  #Forager arguments NA
  scenario$nests[[1]]$whatCurr=NA
  expect_error(optimLoadCurr(u,scenario))

  #Forager arguments missing
  scenario$nests[[1]]$whatCurr=NULL
  expect_error(optimLoadCurr(u,scenario))
})

