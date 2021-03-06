context('Forager moving function (moveForagers) and selection function (whichMoves)')

#Create scenario
nests<-list(xloc=1,yloc=1,n=matrix(c(5,10,15),1),
                       whatCurr="eff",sol=T,
                       eps=0,L_max=59.5,v=7.8,
                       betaVal=0.102/59.5,p_i=1,
                       h=matrix(rep(1.5,3),1),
                       c_f=0.05,c_i=0.0042,
                       H=100,d=matrix(c(50,100,150),1),
                       L=matrix(rep(59.5,3),1),
                       curr=matrix(c(5,10,15),1),
                       steps=5,stepNum=1)
world<-list(mu=matrix(rep(8.33e-05,3),1),
            flDens=matrix(rep(520,3),1),
            e=matrix(rep(14.3,3),1),
            l=matrix(rep(1,3),1),
            f=matrix(rep(0.86,3),1),
            cellSize=10,
            S=matrix(rep(1,3),1),
            alphaVal=matrix(rep(0.013,3),1),
            forageType='random')
baseScen <- list(nests=nests,world=world) #Baseline scenario
use <- c(1,2,3) #Cells to use
temp <- lapply(use,optimLoadCurr,scenario=baseScen) #Get optimL, currency, and S values from baseline
baseScen$world$S[use] <- sapply(temp,function(x) x$S) #Assigns S-values to baseline scenario

for(u in 1:length(temp)){ #For each cell processed
  baseScen$nests$L[use[u]]=temp[[u]][['optimL']][[1]] #Assigns L
  baseScen$nests$curr[use[u]]=temp[[u]][['optimCurr']][[1]] #Assigns currency
}

#Create scenario set
bestScen=makeBest(baseScen)
worstScen=makeWorst(baseScen)
scenSet=list(best=bestScen,base=baseScen,worst=worstScen)

test_that('Currency calculations work properly',{
  #Starting scenario set
  expect_equal(scenSet$base$world$S,matrix(c(0.1409129,0.07139869,0.03995433),1),tol=1e-4) #S
  expect_equal(scenSet$base$nests$n,matrix(c(5,10,15),1)) #n
  expect_equal(scenSet$base$nests$L,matrix(c(0.9601998,0.5816454,0.4577625),1),tol=1e-4) #L
  expect_equal(scenSet$base$nests$curr,matrix(c(9.064411,3.018284,1.277003),1),tol=1e-4) #curr

  expect_equal(scenSet$best$nests$n,matrix(c(10,15,20),1),tol=1e-4) #n in Best scenario
  expect_equal(scenSet$best$nests$curr,matrix(c(4.476088,1.750862,0.7264271),1),tol=1e-4) #curr in Best scenario

  #Solitary foraging case:
  moves=whichMoves(scenSet) #What moves should foragers make?
  expect_equal(list(move=T,from=matrix(c(F,F,T),1),to=matrix(c(T,F,F),1)),moves) #From cell 3 to 1
  newScenSet=moveForagers(scenSet,moves) #Move foragers and save scenario set

  #Since foragers were moved from cell 3 to 1, cell 2 should be identical in both scenarios
  expect_true(newScenSet$base$world$S[2]==scenSet$base$world$S[2]) #S
  expect_true(newScenSet$base$nests$n[2]==scenSet$base$nests$n[2]) #n
  expect_true(newScenSet$base$nests$L[2]==scenSet$base$nests$L[2]) #L
  expect_true(newScenSet$base$nests$curr[2]==scenSet$base$nests$curr[2]) #curr

  #Cells 1 should NOT be equal, since foragers were moved
  expect_false(newScenSet$base$world$S[1]==scenSet$base$world$S[1]) #S
  expect_false(newScenSet$base$nests$n[1]==scenSet$base$nests$n[1]) #n
  expect_false(newScenSet$base$nests$L[1]==scenSet$base$nests$L[1]) #L
  expect_false(newScenSet$base$nests$curr[1]==scenSet$base$nests$curr[1]) #curr

  #Cells 3 should NOT be equal, since foragers were moved
  expect_false(newScenSet$base$world$S[3]==scenSet$base$world$S[3]) #S
  expect_false(newScenSet$base$nests$n[3]==scenSet$base$nests$n[3]) #n
  expect_false(newScenSet$base$nests$L[3]==scenSet$base$nests$L[3]) #L
  expect_false(newScenSet$base$nests$curr[3]==scenSet$base$nests$curr[3]) #curr

  #Move to fixation.
  moves <- whichMoves(newScenSet)
  newScenSet <- moveForagers(newScenSet,moves)
  moves <- whichMoves(newScenSet)

  expect_equal(list(move=F,from=NA,to=NA),moves) #No move. Has reached fixation

  expect_equal(newScenSet$base$world$S,matrix(c(0.08609854,0.07139869,0.068875),1),tol=1e-4)
  expect_equal(newScenSet$base$nests$n,matrix(c(15,10,5),1),tol=1e-4) #Foragers moved to nearest cell
  expect_equal(newScenSet$base$nests$L,matrix(c(0.3194372,0.5816454,1.552095),1),tol=1e-4)
  expect_equal(newScenSet$base$nests$curr,matrix(c(2.793796,3.018284,5.433004),1),tol=1e-4)

})
