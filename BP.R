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
  
  ###################################### Bandt-Pompe functions #############################################
  
  bandt.pompe <- function(series, dimension, delay){
    dyn.load("BandtPompe.so")
    elements = formationPattern(series, dimension, delay, 1)
    element.size = dim(elements)[1]
    probability <- .Call("BandtPompe", elements, dimension, element.size)
    return(probability)
  }
  
  hc.analysis <- function(dimension){
    j = 1
    r = seq(from = 3.5, to = 4.0, by = 0.01)
    Entropy.Complexity = matrix(nrow = length(r), ncol = 2)
    
    for(rr in r){
        series =  series.generator.map(rr)
        g = bandt.pompe(series, dimension, dimension)
        Entropy.Complexity[j, 1] <- shannonNormalized(as.vector(g))
        Entropy.Complexity[j, 2] <- Ccomplexity(as.vector(g))
        j = j + 1
    }
    cat("D: ", dimension, " - tal: ", 1, "\n")
    write.csv(Entropy.Complexity, paste0("Data/BP_HC_D", dimension, "T", dimension, ".csv"))
  }
  
  
  registerDoParallel(cores = 4)
  foreach(i = 3:6) %dopar% {
    hc.analysis(i)
  }
  
  
