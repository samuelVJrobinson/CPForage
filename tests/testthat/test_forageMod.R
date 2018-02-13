context('Full function (forageMod)')

#Create test world for run
nu_i<-0.15/3600 #Nectar production/hr for a single flower
flDens<-520 #Flower density/m2
e_i<-14.35 #Energetic value/unit
l_i<-0.87 #Canola standing crop (0.87uL)
f_i<-0.86 #Inter-flower flight time
#World structure
cellSize<-10 #10m cells (100m^2)
worldSize<-120 #120x120m field (100x100m field with 10m buffer zone worth nothing)
world1<-list(mu=matrix(0,worldSize/cellSize,worldSize/cellSize),  #Empty world
           flDens=matrix(0,worldSize/cellSize,worldSize/cellSize),
           e=matrix(0,worldSize/cellSize,worldSize/cellSize),
           l=matrix(0,worldSize/cellSize,worldSize/cellSize),
           f=matrix(0,worldSize/cellSize,worldSize/cellSize),
           cellSize=cellSize,
           forageType='random') #Competition for flowers within patch
world1$mu[c(2:11),c(2:11)]<-nu_i #Per-flower nectar production in
# canola-filled cells
world1$flDens[c(2:11),c(2:11)]<-flDens #Flower number per cell
world1$e[c(2:11),c(2:11)]<-e_i #Energy production in canola-filled cells
world1$l[c(2:11),c(2:11)]<-l_i #Standing crop in cells with no competition
world1$f[c(2:11),c(2:11)]<-f_i #Inter-flower flight time world1$patchLev=F

#Constants for foragers
honeybeeConstants<-list(L_max=59.5, #Max load capacity (uL)
                      v=7.8, #Velocity (m/s) - Unloaded flight speed
                      beta=0.102, #Proportion reduction in completely loaded flight speed (1-v/v_l)
                      p_i=1, # Max loading rate (uL/s)
                      h=1.5, #Handling time per flower (s)
                      c_f=0.05, #Unloaded flight energetic cost (J/s)
                      c_i=0.0042, #Cost of non-flying activity
                      H=100) #Time spent in the hive (s)

#Nest structure (social rate maximizers)
nests1<-list(nest1=list(xloc=1,yloc=1,n=1000,whatCurr='rat',sol=F,constants=honeybeeConstants,eps=1e-10))

#Run full model (serial)
testOutput1<-forageMod(world1,nests1,2000,verbose=F,parallel=F)

test_that("Results in correct format",{
  expect_length(testOutput1,2) #Nest and world list
  expect_length(testOutput1$nests,1) #Single nest
  expect_length(testOutput1$nests[[1]],22) #22 elements within nest structure
  expect_equal(dim(testOutput1$nests[[1]]$n),c(12,12)) #12 x 12 matrix
})

test_that("Results are consistent",{
  expect_equal(testOutput1$world$S[5,5],0.4253825,tol=1e-4) #S-value
  expect_equal(testOutput1$nests[[1]]$n[5,5],10) #n
  expect_equal(testOutput1$nests[[1]]$L[5,5],59.49993,tol=1e-4) #L
})

test_that('forageMod error handling works',{
  nests2 <- nests1
  nests2[[1]]$xloc <- NULL #Get rid of xloc argument
  expect_error(forageMod(world1,nests2,2000,verbose=F,parallel=F))
  nests2[[1]]$xloc <- 1
  nests2[[1]]$whatCurr <- 'wrong' #Bad currency input
  expect_error(forageMod(world1,nests2,2000,verbose=F,parallel=F))

  world2 <- world1
  world2$mu <- NULL #Get rid of world argument
  expect_error(forageMod(world2,nests1,2000,verbose=F,parallel=F))
})
