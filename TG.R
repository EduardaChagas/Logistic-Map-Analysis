########################################################################################################
# Author: Eduarda Chagas
# Date : Jan 11, 2021
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

######################################Bandt-Pompe functions #############################################

pattern.wedding <- function(patterns){
  #if(max(patterns) == 2)
  #  patterns = patterns + 1
  m = dim(patterns)[1]
  D = dim(patterns)[2]
  symbols = define.symbols(D)
  wedding = rep(0, m)
  for(i in 1:m){
    e = 0
    j = 1
    stop = F
    while(j <= factorial(D) && stop == F){
      if(sum(symbols[j,] == patterns[i,]) == D){
        wedding[i] = j
        stop = T
      }
      j = j + 1
    }
  }
  return(wedding)
}

transition.graph <- function(series, dimension, delay){
  
  graph = matrix(0, nrow = factorial(dimension), ncol = factorial(dimension))
  patterns = formationPattern(series, dimension, delay, 0)
  wedding = pattern.wedding(patterns)
  m = length(wedding)
  
  for(i in 1:(m-1)){
    graph[wedding[i], wedding[i+1]] = graph[wedding[i], wedding[i+1]] + 1
  }
  
  graph = graph/(m-1)
  graph = as.vector(graph)
  return(graph)
}

TG <- function(series, dimension, delay){
  patterns = formationPattern(series, dimension, delay, 0)
  wedding = pattern.wedding(patterns)
  size = length(wedding)
  
  dyn.load("TransitionGraph.so")
  probability <- .Call("TransitionGraph", wedding, dimension, size)
  
  return(probability)
}

hc.analysis <- function(dimension){
  j = 1
  r = seq(from = 3.5, to = 4.0, by = 0.01)
  Entropy.Complexity = matrix(nrow = length(r), ncol = 2)
  
  for(rr in r){
    series =  series.generator.map(rr)
    g = TG(series, dimension, dimension)
    Entropy.Complexity[j, 1] <- shannonNormalized(as.vector(g))
    Entropy.Complexity[j, 2] <- Ccomplexity(as.vector(g))
    j = j + 1
  }
  cat("D: ", dimension, " - tal: ", 1, "\n")
  write.csv(Entropy.Complexity, paste0("Data/TG_HC_D", dimension, "T", dimension, ".csv"))
}


registerDoParallel(cores = 4)
foreach(i = 3:6) %dopar% {
  hc.analysis(i)
}


