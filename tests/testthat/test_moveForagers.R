context('Forager moving function (moveForagers) and selection function (whichMoves)')

#Create scenario
nests<-list(nest1=list(xloc=1,yloc=1,n=matrix(c(5,10,15),1),
                       whatCurr="eff",sol=T,
                       eps=0,L_max=59.5,v=7.8,
                       beta=0.102,p_i=1,
                       h=matrix(rep(1.5,3)),
                       c_f=0.05,c_i=0.0042,
                       H=100,d=matrix(c(50,100,150),1),
                       L=matrix(rep(59.5,3),1),
                       curr=matrix(c(5,10,15),1),
                       steps=5,stepNum=1))
world<-list(mu=matrix(rep(8.33e-05,3),1),
            flDens=matrix(rep(520,3),1),
            e=matrix(rep(14.3,3),1),
            l=matrix(rep(1,3),1),
            f=matrix(rep(0.86,3),1),
            cellSize=10,patchLev=FALSE,
            S=matrix(rep(1,3),1))
baseScen=list(nests=nests,world=world)
use<-c(1,2,3)
temp=lapply(use,optimLoadCurr,scenario=baseScen)
baseScen$world$S[use]=sapply(temp,function(x) x$S) #Assigns S-values
for(u in 1:length(temp)){ #For each cell processed
  for(name in names(temp[[u]]$optimCurr)){ #For each nest within temp[[u]]
  baseScen$nests[[name]][['L']][use[u]]=temp[[u]][['optimL']][[name]] #Assigns L
    baseScen$nests[[name]][['curr']][use[u]]=temp[[u]][['optimCurr']][[name]] #Assigns currency
  }
}
#Create scenario set
bestScen=makeBest(baseScen,1)
worstScen=makeWorst(baseScen,1)
scenSet=list(best=bestScen,base=baseScen,worst=worstScen)

test_that('Currency calculations work properly',{
  #Starting scenario set
  expect_equal(scenSet$base$world$S,matrix(c(0.8170799,0.668982,0.5572836),1),tol=1e-4) #S
  expect_equal(scenSet$base$nests$nest1$n,matrix(c(5,10,15),1)) #n
  expect_equal(scenSet$base$nests$nest1$L,matrix(c(46.81022,51.86594,54.28131),1),tol=1e-4) #L
  expect_equal(scenSet$base$nests$nest1$curr,matrix(c(237.3381,169.9,131.9003),1),tol=1e-4) #curr

  expect_equal(scenSet$best$nests$nest1$n,matrix(c(10,15,20),1),tol=1e-4) #n in Best scenario
  expect_equal(scenSet$best$nests$nest1$curr,matrix(c(216.3022,156.0375,122.1375),1),tol=1e-4) #n in Best scenario

  #Solitary foraging case:
  moves=whichMoves(scenSet,1) #What moves should foragers make?
  # expect_equal(list(move=T,from=matrix(c(F,F,T),1),to=matrix(c(T,F,F),1)),moves) #From cell 3 to 1
  newScenSet=moveForagers(scenSet,1,moves) #Move foragers and save scenario set

  #Since foragers were moved from cell 3 to 1, cell 2 should be identical in both scenarios
  expect_true(newScenSet$base$world$S[2]==scenSet$base$world$S[2]) #S
  expect_true(newScenSet$base$nests$nest1$n[2]==scenSet$base$nests$nest1$n[2]) #n
  expect_true(newScenSet$base$nests$nest1$L[2]==scenSet$base$nests$nest1$L[2]) #L
  expect_true(newScenSet$base$nests$nest1$curr[2]==scenSet$base$nests$nest1$curr[2]) #curr

  #Cells 1 should NOT be equal, since foragers were moved
  expect_false(newScenSet$base$world$S[1]==scenSet$base$world$S[1]) #S
  expect_false(newScenSet$base$nests$nest1$n[1]==scenSet$base$nests$nest1$n[1]) #n
  expect_false(newScenSet$base$nests$nest1$L[1]==scenSet$base$nests$nest1$L[1]) #L
  expect_false(newScenSet$base$nests$nest1$curr[1]==scenSet$base$nests$nest1$curr[1]) #curr

  #Cells 3 should NOT be equal, since foragers were moved
  expect_false(newScenSet$base$world$S[3]==scenSet$base$world$S[3]) #S
  expect_false(newScenSet$base$nests$nest1$n[3]==scenSet$base$nests$nest1$n[3]) #n
  expect_false(newScenSet$base$nests$nest1$L[3]==scenSet$base$nests$nest1$L[3]) #L
  expect_false(newScenSet$base$nests$nest1$curr[3]==scenSet$base$nests$nest1$curr[3]) #curr

  #Move to fixation.
  moves=whichMoves(newScenSet,1) #What moves should foragers make?
  expect_equal(list(move=T,from=matrix(c(F,F,T),1),to=matrix(c(T,F,F),1)),moves) #From cell 3 to 1
  newScenSet=moveForagers(newScenSet,1,moves) #Move foragers and save scenario set
  moves=whichMoves(newScenSet,1) #What moves should foragers make?
  expect_equal(list(move=T,from=matrix(c(F,F,T),1),to=matrix(c(T,F,F),1)),moves) #From cell 3 to 1
  newScenSet=moveForagers(newScenSet,1,moves) #Move foragers and save scenario set
  moves=whichMoves(newScenSet,1) #What moves should foragers make?
  expect_equal(list(move=F,from=NA,to=NA),moves) #Done

})
