########################################################################################################
# Author: Eduarda Chagas
# Date : Jan 21, 2021
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################

############################################# Packages #################################################
if(!require(gtools)) install.packages("gtools")
if(!require(doParallel)){
  install.packages("doParallel")
  require(doParallel)
}
source("theory_information.R")
source("BandtPompeAux.R")
source("LogisticMap.R")

###################################### AAPE functions #############################################

AAPE <- function(series, dimension, delay, A = 0.5){
  series = unlist(series)
  symbols = formationPattern(series,dimension,delay,0)
  elements = formationPattern(series,dimension,delay,1)
  patterns = define.symbols(dimension)
  sw = rep(0,factorial(dimension))
  
  for(i in 1:factorial(dimension)){
    for(j in 1:dim(symbols)[1]){
      if(all(symbols[j,] != patterns[i,])){
        sw[i] = sw[i] + ((A/dimension)*abs(elements[j,1])) + sum((A/dimension)*abs(elements[j,2:dimension ]) + ((1-A)/(dimension-1))*abs(elements[j,2:dimension] - elements[j,1:(dimension - 1)]))
      }
    }
  }
  
  pw = sw/sum(sw)
  pw
}

aape.analysis <- function(dimension){
  j = 1
  r = seq(from = 3.5, to = 4.0, by = 0.01)
  Entropy.Complexity = matrix(nrow = length(r), ncol = 2)
  
  for(rr in r){
    series =  series.generator.map(rr)
    g = AAPE(series, dimension, 1)
    Entropy.Complexity[j, 1] = shannonNormalized(as.vector(g))
    Entropy.Complexity[j, 2] = Ccomplexity(as.vector(g))
    j = j + 1
  }
  cat("D: ", dimension, " - tal: ", 1, "\n")
  write.csv(Entropy.Complexity, paste0("Data/AAPE_HC_D", dimension, "T", 1, ".csv"))
}


registerDoParallel(cores = 4)
foreach(i = 3:6) %dopar% {
  aape.analysis(i)
}


