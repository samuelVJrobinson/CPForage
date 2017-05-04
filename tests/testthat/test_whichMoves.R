context('Movement selection function (whichMoves)')

#Reduced nest structure with n, curr, sol, stepNum, and steps (tranfer number)
nests<-list(nest1=list(n=matrix(c(0,5,10),nrow=1),curr=matrix(c(11,10,9),nrow=1),sol=T,stepNum=1,steps=5))
world<-NA #Empty world
base<-list(nests=nests,world=world) #Base scenario
scenSet<-list(best=base,
              base=base,
              worst=base) #Scenario set
scenSet$best$nests$nest1$n=scenSet$best$nests$nest1$n+5 #Number in "Best nests"
scenSet$best$nests$nest1$curr=scenSet$best$nests$nest1$curr-1 #Curr in "Best nests"
scenSet$worst$nests$nest1$n=matrix(c(0,0,5),nrow=1) #Number in "worst nests"
scenSet$worst$nests$nest1$curr=matrix(c(11,11,10),nrow=1) #Curr in "worst nests"

test_that('Move calculations work properly',{
  moves<-whichMoves(scenSet,i=1) #Best move for solitary foragers
  expect_equal(which(moves$from),3)
  expect_equal(which(moves$to),1)

  scenSet$best$nests$nest1$sol=F
  scenSet$base$nests$nest1$sol=F
  scenSet$worst$nests$nest1$sol=F
  moves<-whichMoves(scenSet,i=1) #Best move for social foragers
  expect_equal(which(moves$from),3)
  expect_equal(which(moves$to),1)

  scenSet$best$nests$nest1$curr[1]=3
  moves<-whichMoves(scenSet,i=1)
  expect_equal(which(moves$to),2)
})




