context('Summed currency selection (curr)')

params <- list(L_i=59.5,
            L_max_i=59.5,
            n_i=10, #10 foragers
            h_i=1.5,
            p_i=1,
            f_i=0.86,
            d_i=100, #100m
            v_i=7.8,
            beta_i=0.102/59.5,
            H_i=100,
            c_f=0.05,
            c_i=0.0042,
            mu=0.3/3600,
            l=1,
            e=14.35,
            NumFls=520, #520 flowers
            alphaVal=0.013)

# #Plot of function results
# par(mfrow=c(3,1))
# N <- seq(1,100000,10)
# temp <- t(sapply(N,function(x){
#   with(params,curr(L_i,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                    c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,
#                    alphaVal=alphaVal,forageType='omniscient'))
# }))
# plot(N,temp[,1],type='l',xlab='N foragers',ylab='Eff')
# temp <- t(sapply(N,function(x){
#   with(params,curr(L_i,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                    c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,
#                    alphaVal=alphaVal,forageType='omniscient'))
# }))
# plot(N,temp[,1],type='l',xlab='N foragers',ylab='Rat')
# plot(N,temp[,2],type='l',xlab='N foragers',ylab='S')

test_that('Summed currency works properly',{
  #Efficiency, using omniscient foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=T,
                                alphaVal=alphaVal,forageType='omniscient')),
               1.174066,tol=1e-4)
  #Efficiency, using random foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                c_i,c_f,mu,l,e,NumFls,alphaVal=alphaVal,
                                whatCurr_i='eff',sumAll=T,forageType='random')),
               1.151189,tol=1e-4)
  #Efficiency, using NN foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                c_i,c_f,mu,l,e,NumFls,alphaVal=alphaVal,
                                whatCurr_i='eff',sumAll=T,forageType='nn')),
               1.151189,tol=1e-4)

  #Net rate, using omniscient foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                c_i,c_f,mu,l,e,NumFls,alphaVal=alphaVal,
                                whatCurr_i='rat',sumAll=T,forageType='omniscient')),
               0.03358101,tol=1e-4)
  #Net rate, using random foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                c_i,c_f,mu,l,e,NumFls,alphaVal=alphaVal,
                                whatCurr_i='rat',sumAll=T,forageType='random')),
               0.03292918,tol=1e-4)
  #Net rate, using NN foraging
  expect_equal(with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                                c_i,c_f,mu,l,e,NumFls,alphaVal=alphaVal,
                                whatCurr_i='rat',sumAll=T,forageType='nn')),
               0.03292918,tol=1e-4)
})

test_that('Individual currency and S-values work properly',{
  #Efficiency with 1 nest
  effRes1 <- with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                        c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,forageType='random',
                        alphaVal=alphaVal))
  expect_equal(length(effRes1),2)
  expect_equal(effRes1[['S']],0.0104476)
  expect_equal(effRes1[['eff']],1.151189)

  #Same thing, but 200 foragers
  effRes2 <- with(params,curr(L_i,L_max_i,n_i=200,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,c_i,c_f,
                        mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,forageType='random',alphaVal=alphaVal))
  expect_equal(length(effRes2),2)
  expect_equal(effRes2[['S']],0.0005209168)
  expect_equal(effRes2[['eff']],-0.892096,tol=1e-06)

  #Net rate with 1 nests
  ratRes1 <- with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                        c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='random',alphaVal=alphaVal))
  expect_equal(length(ratRes1),2)
  expect_equal(ratRes1[['S']],0.0104476)
  expect_equal(ratRes1[['rat']],0.03294044)

  #Same thing, but 200 foragers
  ratRes2 <- with(params,curr(L_i,L_max_i,200,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                        c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='random',alphaVal=alphaVal))
  expect_equal(length(ratRes2),2)
  expect_equal(ratRes2[['S']],0.0005209168)
  expect_equal(ratRes2[['rat']],-0.02569135)

  #NN foraging:
  #Net rate with 1 nests
  ratRes1 <- with(params,curr(L_i,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                        c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='nn',alphaVal=alphaVal))
  expect_equal(length(ratRes1),2)
  expect_equal(ratRes1[['S']],0.0104476)
  expect_equal(ratRes1[['rat']],0.03294044)

  #Same thing, but 200 foragers
  ratRes2 <- with(params,curr(L_i,L_max_i,200,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                              c_i,c_f,mu,l,e,NumFls,whatCurr_i='rat',sumAll=F,forageType='nn',alphaVal=alphaVal))
  expect_equal(length(ratRes2),2)
  expect_equal(ratRes2[['S']],0.0005209168)
  expect_equal(ratRes2[['rat']],-0.02569135)

  #NN foraging, but with large numbers of flowers (90000) relative to foragers (3). S should be 1.
  ratRes3 <- with(params,curr(L_i,L_max_i,n_i=3,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                        c_i,c_f,mu,l,e,NumFls=400*15^2,whatCurr_i='rat',sumAll=F,forageType='nn',alphaVal=alphaVal))
  expect_equal(length(ratRes3),2)
  expect_equal(ratRes3[['S']],1)
  expect_equal(ratRes3[['rat']],2.5889,tol=1e-5)

  #However, this should not be true with random foragers. S should be 0.92.
  ratRes4 <- with(params,curr(L_i,L_max_i,n_i=3,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                  c_i,c_f,mu,l,e,NumFls=400*15^2,whatCurr_i='rat',sumAll=F,forageType='random',alphaVal=alphaVal))
  expect_equal(length(ratRes4),2)
  expect_equal(ratRes4[['S']],0.9298684)
  expect_equal(ratRes4[['rat']],2.505374,tol=1e-5)
})


test_that('Exception handling (limited - mostly in optimLoadCurr)',{
  #Load of 0
  effRes1 <- with(params,curr(0,L_max_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
                        c_i,c_f,mu,l,e,NumFls,whatCurr_i='eff',sumAll=F,forageType='random',alphaVal=alphaVal))
  expect_equal(effRes1[['S']],1)
  expect_equal(effRes1[['eff']],-1)
})

# #Demo of NN foraging with n flowers vs random foraging
# par(mfrow=c(2,1))
# for(ft in c('nn','random')){
#   plot(seq(1000,90000,1000),sapply(seq(1000,90000,1000),function(x) {
#     with(params,curr(L_i,L_max_i,n_i=4,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                      c_i,c_f,mu,l,e,NumFls=x,whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#   }),xlab='n flowers',ylab='S',type='p',pch=19,main=paste(ft,'foraging',sep=' '),ylim=c(0,1))
#   points(seq(1000,90000,1000),sapply(seq(1000,90000,1000),function(x) {
#     with(params,curr(L_i,L_max_i,n_i=3,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                      c_i,c_f,mu,l,e,NumFls=x,whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#   }),col='blue',pch=19)
#   points(seq(1000,90000,1000),sapply(seq(1000,90000,1000),function(x) {
#     with(params,curr(L_i,L_max_i,n_i=2,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                      c_i,c_f,mu,l,e,NumFls=x,whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#   }),col='green',pch=19)
#   abline(h=1,lty='dashed')
# }
# legend('bottomright',c('4','3','2'),fill=c('black','blue','green'),title='Number of foragers')
#
# #Demo of NN foraging with n foragers, across varying load sizes and flower numbers
# par(mfrow=c(2,3))
# for(ft in c('nn','random')){
#   for(div in 1:3){
#     plot(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=400*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),xlab='n foragers',ylab='S',type='p',pch=19,ylim=c(0,1),
#     main=paste(ft,'foraging,','load/',div,sep=' '))
#     points(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=200*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),col='blue',pch=19)
#     points(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=50*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),col='green',pch=19)
#     points(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=5*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),col='red',pch=19)
#
#     abline(h=1,lty='dashed')
#   }
# }
# legend('topright',c('400','200','50','5'),fill=c('black','blue','green','red'),title='Flowers per m2')
#
# #Same thing, but with higher mu and l
#
# params$mu <- 1/3600 #1 uL/hr
# params$l <- 3 #3 uL max nectar
#
# par(mfrow=c(2,1))
# for(ft in c('nn','random')){
#   plot(seq(1000,90000,1000),sapply(seq(1000,90000,1000),function(x) {
#     with(params,curr(L_i,L_max_i,n_i=4,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                      c_i,c_f,mu,l,e,NumFls=x,whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#   }),xlab='n flowers',ylab='S',type='p',pch=19,main=paste(ft,'foraging',sep=' '),ylim=c(0,1))
#   points(seq(1000,90000,1000),sapply(seq(1000,90000,1000),function(x) {
#     with(params,curr(L_i,L_max_i,n_i=3,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                      c_i,c_f,mu,l,e,NumFls=x,whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#   }),col='blue',pch=19)
#   points(seq(1000,90000,1000),sapply(seq(1000,90000,1000),function(x) {
#     with(params,curr(L_i,L_max_i,n_i=2,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                      c_i,c_f,mu,l,e,NumFls=x,whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#   }),col='green',pch=19)
#   abline(h=1,lty='dashed')
# }
# legend('bottomright',c('4','3','2'),fill=c('black','blue','green'),title='Number of foragers')
#
#
# par(mfrow=c(2,3))
# for(ft in c('nn','random')){
#   for(div in 1:3){
#     plot(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=400*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),xlab='n foragers',ylab='S',type='p',pch=19,ylim=c(0,1),
#     main=paste(ft,'foraging,','load/',div,sep=' '))
#     points(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=200*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),col='blue',pch=19)
#     points(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=50*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),col='green',pch=19)
#     points(seq(1,100,1),sapply(seq(1,100,1),function(x) {
#       with(params,curr(L_i/div,L_max_i,n_i=x,h_i,p_i,f_i,d_i,v_i,beta_i,H_i,
#                        c_i,c_f,mu,l,e,NumFls=5*(10^2),whatCurr_i='rat',sumAll=F,forageType=ft))[2]
#     }),col='red',pch=19)
#
#     abline(h=1,lty='dashed')
#   }
# }
# legend('topright',c('400','200','50','5'),fill=c('black','blue','green','red'),title='Flowers per m2')

