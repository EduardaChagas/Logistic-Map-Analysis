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

weights <- function(series, dimension, delay){
  groups = formationPattern(series,dimension,delay,1)
  weight = w = rep(0, dim(groups)[1])
  for(i in 1:dim(groups)[1]){
    weight[i] = (sum((groups[i,] - mean(groups[i,]))^2))/dimension
  }
  weight
}

WPE <- function(series, dimension, delay){
  series = unlist(series)
  w = weights(series,dimension,delay)
  symbols = formationPattern(series,dimension,delay,0)
  patterns = define.symbols(dimension)
  sw = rep(0,factorial(dimension))
  
  for(i in 1:factorial(dimension)){
    for(j in 1:dim(symbols)[1]){
      if(all(symbols[j,] != patterns[i,]))
        sw[i] = sw[i] + w[j]
    }
  }
  pw = sw/sum(sw)
  pw
}

wpe.analysis <- function(dimension){
  j = 1
  r = seq(from = 3.5, to = 4.0, by = 0.01)
  Entropy.Complexity = matrix(nrow = length(r), ncol = 2)
  
  for(rr in r){
    series =  series.generator.map(rr)
    g = WPE(series, dimension, 1)
    Entropy.Complexity[j, 1] = shannonNormalized(as.vector(g))
    Entropy.Complexity[j, 2] = Ccomplexity(as.vector(g))
    j = j + 1
  }
  cat("D: ", dimension, " - tal: ", 1, "\n")
  write.csv(Entropy.Complexity, paste0("Data/WPE_HC_D", dimension, "T", 1, ".csv"))
}


registerDoParallel(cores = 4)
foreach(i = 3:6) %dopar% {
  wpe.analysis(i)
}


