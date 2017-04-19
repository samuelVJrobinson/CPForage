#Convert model runs to dataframes
nests2df=function(nests){
  require(dplyr)
  require(tidyr)
  convertNest=function(i,nes){
    nNests=bind_rows(data.frame(nes[[i]]$n,rows=0:(nrow(nes[[i]]$n)-1),category='numForagers'), #Occupancy
                     data.frame(nes[[i]]$L,rows=0:(nrow(nes[[i]]$n)-1),category='loadSize'), #Load size
                     data.frame(nes[[i]]$curr,rows=0:(nrow(nes[[i]]$n)-1),category='currency'), #Currency
                     data.frame(nes[[i]]$d,rows=0:(nrow(nes[[i]]$n)-1),category='distance'), #Distance to patch
                     data.frame(nes[[i]]$loadingTime,rows=0:(nrow(nes[[i]]$n)-1),category='loadingTime'), #Time spent loading
                     data.frame(nes[[i]]$travelTime,rows=0:(nrow(nes[[i]]$n)-1),category='travelTime'), #Travel time to patch
                     data.frame(nes[[i]]$boutLength,rows=0:(nrow(nes[[i]]$n)-1),category='boutLength')) #Total length of foraging bout
    names(nNests)[1:nrow(nes[[i]]$n)]=c(0:(nrow(nes[[i]]$n)-1))
    nNests=gather(nNests,cols,value,-rows:-category,convert=T) %>%
      spread(category,value) %>%
      mutate(nest=i,rate=nes[[i]]$rate,sol=nes[[i]]$sol)
    return(nNests)
  }
  return(bind_rows(lapply(1:length(nests),nes=nests,convertNest)))
}
