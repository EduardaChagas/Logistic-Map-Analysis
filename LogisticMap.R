########################################################################################################
# Author: Eduarda Chagas
# Date : Jan 23, 2021
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################


series.generator.map <- function(r){
  time.series <- vector(mode="numeric", length = 10^6)
  index = 1
  time.series[index] <- 0.1
  for(i in 2:(10000 + (10^6))){
    time.series[index] <- r * time.series[index] * (1 - time.series[index])
    if(i > 10001){
      index = index + 1
      time.series[index] <- r * time.series[index-1] * (1 - time.series[index-1])
    }
  }
  return(time.series)
}