#Convert model runs to dataframes
nests2df <- function(modelRun){
  require(dplyr)
  require(tidyr)
  convertNest=function(i,nes){
    nNests=bind_rows(data.frame(nes[[i]]$n,rows=0:(nrow(nes[[i]]$n)-1),category='numForagers'), #Occupancy
               data.frame(nes[[i]]$L,rows=0:(nrow(nes[[i]]$n)-1),category='loadSize'), #Load size
               data.frame(nes[[i]]$curr,rows=0:(nrow(nes[[i]]$n)-1),category='currency'), #Currency
               data.frame(nes[[i]]$d,rows=0:(nrow(nes[[i]]$n)-1),category='distance'), #Distance to patch
               data.frame(nes[[i]]$loadingTime,rows=0:(nrow(nes[[i]]$n)-1),category='loadingTime'), #Time spent loading
               data.frame(nes[[i]]$travelTime,rows=0:(nrow(nes[[i]]$n)-1),category='travelTime'), #Travel time to patch
               data.frame(nes[[i]]$boutLength,rows=0:(nrow(nes[[i]]$n)-1),category='boutLength')) %>% #Length of foraging bout
    setNames(gsub('X','',names(.))) %>% #Remove "X" from names
    gather(cols,value,-rows:-category,convert=T) %>% #Gather into single column
    spread(category,value) %>% #Spread over categories
    mutate(nest=i,whichCurr=nes[[i]]$whatCurr,sol=nes[[i]]$sol) #Add nest number, currency, and solitary foraging
    return(nNests)
  }

  #Warnings about factor to character conversion suppressed
  nests=suppressWarnings(bind_rows(lapply(1:length(modelRun$nests),nes=modelRun$nests,convertNest)))
  world=data.frame(modelRun$world$S,rows=0:(nrow(modelRun$world$S)-1),category='S') %>%
    setNames(gsub('X','',names(.))) %>%
    gather(cols,S,-rows:-category) %>%
    select(-category) %>%
    mutate(cols=as.numeric(cols)-1)
  return(list(nests=nests,world=world))
}
