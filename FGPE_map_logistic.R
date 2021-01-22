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

###################################### AAPE functions #############################################

define.symbols <- function(dimension){
  d = c(1:dimension)
  symbol = matrix(unlist(permutations(n=dimension,r=dimension,v=d)),nrow = factorial(dimension),ncol = dimension,byrow = FALSE)
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

FGPE <- function(series, dimension, delay, p.alpha = 0.5){
  series = unlist(series)
  patterns = formationPattern(series, dimension, delay, 0)
  elements = formationPattern(series, dimension, delay, 1)
  n.patterns = matrix(nrow = dim(patterns)[1], ncol = dimension + 1)
  
  for(i in 1:dim(patterns)[1]){
    dj = abs(elements[i, 2:dimension] - elements[i, 1:(dimension-1)])
    if(sd(dj) != 0)
      q = floor(max(dj)/(sd(dj)*p.alpha))
    else
      q = 0
    n.patterns[i, 1:dimension] = patterns[i,]
    n.patterns[i, dimension+1] = q
  }
  probability = find.new.pattern(n.patterns, dim(patterns)[1], (dimension+1))
  return(probability)
}

find.new.pattern <- function(n.patterns, n.row, n.col){
  found.patterns = 0
  ci = rep(0, n.row)
  pi.symbols = matrix(nrow = n.row, ncol = n.col)
  
  for(i in 1:n.row){
    found = FALSE
    if(found.patterns == 0){
      found = TRUE
      ci[1] = 1
      found.patterns = 1
      pi.symbols[1,] = n.patterns[i,]
    }
    else{
      for(j in 1:found.patterns){
        if(all(pi.symbols[j,] == n.patterns[i,])){
          found = TRUE 
          ci[j] = ci[j] + 1
          break
        }
      }
    }
    
    if(!found){
      found.patterns = found.patterns + 1
      ci[found.patterns] = 1
      pi.symbols[found.patterns,] = n.patterns[i,]
    } 
  }
  ci = na.omit(ci)
  ci = ci[1:found.patterns]
  pi = ci/sum(ci)
  return(pi)   
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

fgpe.analysis <- function(dimension){
  j = 1
  r = seq(from = 3.5, to = 4.0, by = 0.01)
  Entropy.Complexity = matrix(nrow = length(r), ncol = 2)
  
  for(rr in r){
    series =  series.generator.map(rr)
    g = FGPE(series, dimension, 1)
    Entropy.Complexity[j, 1] = shannonNormalized(as.vector(g))
    Entropy.Complexity[j, 2] = Ccomplexity(as.vector(g))
    j = j + 1
  }
  cat("D: ", dimension, " - tal: ", 1, "\n")
  write.csv(Entropy.Complexity, paste0("Data/FGPE_HC_D", dimension, "T", 1, ".csv"))
}


registerDoParallel(cores = 4)
foreach(i = 3:6) %dopar% {
  fgpe.analysis(i)
}


