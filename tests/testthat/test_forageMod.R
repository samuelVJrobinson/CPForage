context('Full function (forageMod)')

#Create test world for run
nu_i<-0.15/3600 #Nectar production/hr for a single flower
flDens<-200 #Flower density/m2 - represents early canola bloom
e_i<-14.35 #Energetic value/unit
l_i<-0.87 #Canola standing crop (0.87uL)
f_i<-3 #Inter-flower flight time
#World structure
cellSize<-10 #10m cells (100m^2)
worldSize<-120 #120x120m field (100x100m field with 10m buffer zone with weeds in it)
world1<-list(mu=matrix(1/3600,worldSize/cellSize,worldSize/cellSize),  #Weeds produce 1uL/hr
           flDens=matrix(3*(cellSize^2),worldSize/cellSize,worldSize/cellSize), #1 fls/m2
           e=matrix(8,worldSize/cellSize,worldSize/cellSize), #worth about 75% of canola
           l=matrix(3,worldSize/cellSize,worldSize/cellSize), #max nectar = 3uL
           f=matrix(5,worldSize/cellSize,worldSize/cellSize), #5 second flight time b/w fls
           alphaVal=matrix(0.013,worldSize/cellSize,worldSize/cellSize), #Alpha value (J/s*uL)
           cellSize=cellSize,
           forageType='random') #Competition for flowers within patch
world1$mu[c(2:11),c(2:11)]<-nu_i #Per-flower nectar production in
# canola-filled cells
world1$flDens[c(2:11),c(2:11)]<-flDens*cellSize^2 #Flower number per cell
world1$e[c(2:11),c(2:11)]<-e_i #Energy production in canola-filled cells
world1$l[c(2:11),c(2:11)]<-l_i #Standing crop in cells with no competition
world1$f[c(2:11),c(2:11)]<-f_i #Inter-flower flight time world1$patchLev=F

#Constants for foragers
honeybeeConstants<-list(L_max=59.5, #Max load capacity (uL)
                      v=7.8, #Velocity (m/s) - Unloaded flight speed
                      betaVal=0.102/59.5, #Proportion reduction in completely loaded flight speed (1-v/v_l)
                      p_i=1, # Max loading rate (uL/s)
                      h=1.5, #Handling time per flower (s)
                      c_f=0.05, #Unloaded flight energetic cost (J/s)
                      c_i=0.0042, #Cost of non-flying activity
                      H=100) #Time spent in the hive (s)


#Nest structure (social rate maximizers)
nests1<-list(xloc=1,yloc=1,n=1000,whatCurr='rat',sol=F,
                        constants=honeybeeConstants,eps=0,steps=c(50,5,1))
#Run model in serial
testOutput1<-forageMod(world1,nests1,2000,verbose=F,parallel=F)

#Nest structure (social efficiency maximizers)
nests2<-list(xloc=1,yloc=1,n=1000,whatCurr='eff',sol=F,
             constants=honeybeeConstants,eps=0,steps=c(50,5,1))

testOutput2<-forageMod(world1,nests2,2000,verbose=F,parallel=F)

# #Plot of results
# par(mfrow=c(3,1))
# plot(diag(testOutput2$nests$d[2:11,2:11]),diag(testOutput2$world$S[2:11,2:11]),xlab='Distance',ylab='S',pch=19,col='red',cex=1.3,main='Depletion',ylim=c(0,1))
# points(diag(testOutput1$nests$d[2:11,2:11]),diag(testOutput1$world$S[2:11,2:11]),pch=19)
# abline(h=c(0,1),lty='dashed')
# legend('topright',c('Efficiency','Net Rate'),fill=c('red','black'))
#
# plot(diag(testOutput2$nests$d[2:11,2:11]),diag(testOutput2$nests$n[2:11,2:11]),xlab='Distance',ylab='Count',pch=19,col='red',cex=1.3,main='Forager number',
#      ylim=c(0,max(c(testOutput2$nests$n[2:11,2:11],testOutput1$nests$n[2:11,2:11]),na.rm=T)))
# points(diag(testOutput1$nests$d[2:11,2:11]),diag(testOutput1$nests$n[2:11,2:11]),pch=19)
#
# plot(diag(testOutput2$nests$d[2:11,2:11]),diag(testOutput2$nests$L[2:11,2:11]),xlab='Distance',ylab='L',pch=19,col='red',cex=1.3,main='Load size',
#      ylim=c(0,max(c(testOutput2$nests$L,testOutput1$nests$L),na.rm=T)))
# points(diag(testOutput1$nests$d[2:11,2:11]),diag(testOutput1$nests$L[2:11,2:11]),pch=19)
#
# #Heatmaps
#
# #NOTE: foragers in netrate simulation "pile up" in nest cell, because total N foragers is very large compared to size of landscape
# par(mfrow=c(4,2))
# image(testOutput1$world$S,main='NetRate: S'); image(testOutput2$world$S,main='Eff: S')
# image(testOutput1$nests$n,main='NetRate: n'); image(testOutput2$nests$n,main='Eff: n')
# image(testOutput1$nests$L,main='NetRate: L'); image(testOutput2$nests$L,main='Eff: L')
# image(testOutput1$nests$curr,main='NetRate: curr'); image(testOutput2$nests$curr,main='Eff: curr')
# par(mfrow=c(1,1))

test_that("Results in expected format",{
  expect_length(testOutput1,2) #Nest and world list
  expect_length(testOutput1$nests,22) #23 elements within nest structure
  expect_equal(dim(testOutput1$nests$n),c(12,12)) #12 x 12 matrix
})

test_that("Results are consistent",{
  #World 1
  expect_equal(testOutput1$world$S[5,5],0.3462,tol=1e-4) #S-value
  expect_equal(testOutput1$nests$n[5,5],10) #n
  expect_equal(testOutput1$nests$L[5,5],59.5,tol=1e-4) #L

  #World 2
  expect_equal(testOutput2$world$S[5,5],0.5309762,tol=1e-4) #S-value
  expect_equal(testOutput2$nests$n[5,5],11) #n
  expect_equal(testOutput2$nests$L[5,5],6.498771,tol=1e-4) #L
})

test_that('forageMod error handling works',{
  nests2 <- nests1
  nests2$xloc <- NULL #Get rid of xloc argument
  expect_error(forageMod(world1,nests2,2000,verbose=F,parallel=F))
  nests2$xloc <- 1
  nests2$whatCurr <- 'wrong' #Bad currency input
  expect_error(forageMod(world1,nests2,2000,verbose=F,parallel=F))

  world2 <- world1
  world2$mu <- NULL #Get rid of world argument
  expect_error(forageMod(world2,nests1,2000,verbose=F,parallel=F))
})
