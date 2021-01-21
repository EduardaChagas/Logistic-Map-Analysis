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

######################################Bandt-Pompe functions #############################################

define.symbols <- function(dimension){
  d = c(1:dimension)
  symbol = matrix(unlist(permutations(n=dimension,r = dimension)),nrow = factorial(dimension),ncol = dimension,byrow = FALSE)
  symbol = symbol - 1
  symbol
}


FP <- function(n, dimension, delay){
  dyn.load("FormationPatterns.so")
  p <- .Call("FormationPatterns", n, dimension, delay)
  p = t(p) + 1
  return(p)
}

formationPattern <- function(serie, dimension, delay, option){
  i = 1
  n = length(serie)
  p_patterns = elements = index2 = matrix(nrow=n,ncol=dimension)
  index = c(0:(dimension-1))
  
  index2 = FP(length(serie), dimension, delay)
  
  while((i + ((dimension-1)*delay)) <= n){ 
    elements[i,] = serie[index2[i,]]
    p_patterns[i,] = index[order(elements[i,])]
    i = i + 1
  }
  
  if(option == 0){
    p_patterns = na.omit(p_patterns)
    return(p_patterns[1:dim(p_patterns)[1],])
  }else if(option == 1){
    elements = na.omit(elements)
    return(elements[1:dim(elements)[1],])    
  }else{
    index2 = na.omit(index2)
    return(index2[1:dim(index2)[1],])    
  }
}

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

hc.analysis <- function(dimension){
  j = 1
  r = seq(from = 3.5, to = 4.0, by = 0.01)
  Entropy.Complexity = matrix(nrow = length(r), ncol = 2)
  
  for(rr in r){
    series =  series.generator.map(rr)
    g = TG(series, dimension, 1)
    Entropy.Complexity[j, 1] <- shannonNormalized(as.vector(g))
    Entropy.Complexity[j, 2] <- Ccomplexity(as.vector(g))
    j = j + 1
  }
  cat("D: ", dimension, " - tal: ", 1, "\n")
  write.csv(Entropy.Complexity, paste0("Data/TG_HC_D", dimension, "T", 1, ".csv"))
}


registerDoParallel(cores = 4)
foreach(i = 3:6) %dopar% {
  hc.analysis(i)
}


