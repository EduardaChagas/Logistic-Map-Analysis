########################################################################################################
# Author: Eduarda Chagas
# Date : Jan 2021
# Contact: eduarda.chagas@dcc.ufmg.br
########################################################################################################

############################################# Packages #################################################

if(!require(gtools)) install.packages("gtools")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggthemes)) install.packages("ggthemes")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(latex2exp)) install.packages("latex2exp")
source("theory_information.R")

###################################### Function of Plot ################################################

HC.Plane.no.cota <- function(dimension, signal.values){
  
  #XMIN = min(signal.values[,1]) + 0.0005
  #XMAX = min(max(signal.values[,1]) + 0.0005, 1)
  #YMIN = max(0,min(signal.values[,2]))
  #YMAX = max(signal.values[,2])
  signal.values = data.frame("H" = signal.values[,1], "C" = signal.values[,2])
  
  p = cotas(dimension)
  p = p + 
    geom_point(data = signal.values, aes(x = H, y = C), size = 2) +
    labs(x = TeX("\\textit{H}"), y = TeX("\\textit{C}"))  +
    scale_shape_identity() +
    #xlim(limits=c(XMIN, XMAX)) + ylim(limits=c(YMIN, YMAX)) + 
    theme_few(base_size = 14, base_family = "serif")  + 
    theme(plot.title = element_text(hjust=0.5)) + 
    scale_colour_few("Dark")
  print(p)
  return(p)
}

plot.TG.analysis <- function(){
  
  plots = array(list(), 4)
  n.total = 51
  
  d = tau = 3 
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D3T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[1]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 3, Delay = 1")))
  
  d = tau = 4
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D4T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[2]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 4, Delay = 1")))
  
  d = tau = 5
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D5T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[3]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 5, Delay = 1")))
  
  d = tau = 6
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D6T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[4]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 6, Delay = 1")))
  
  pdf("HCAnalysis.pdf", width = 10, height = 7) 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + 
    labs(title = "Ordinal Patterns Transition Graphs") +           
    labs(x = TeX("\\textit{Normalized Entropy}"), y = TeX("\\textit{Statistical Complexity}"))  +
    theme_igray() + theme(text = element_text(size = 16, family="Times", face="italic"), plot.title = element_text(hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  dev.off() 
}
