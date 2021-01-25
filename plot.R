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

HC.Plane.cota <- function(dimension, signal.values){
  signal.values = data.frame("H" = signal.values[,1], "C" = signal.values[,2], "Rho" = factor(c(rep("3.5 - 3.99", 50), "4.0")))
  
  p = cotas(factorial(dimension)^2)
  p = p + 
    geom_point(data = signal.values, aes(x = H, y = C, color = Rho), size = 2) +
    labs(x = TeX("\\textit{H}"), y = TeX("\\textit{C}"))  +
    scale_shape_identity() +
    theme_few(base_size = 14, base_family = "serif")  + 
    theme(plot.title = element_text(hjust=0.5)) + 
    scale_colour_few("Dark")
  print(p)
  return(p)
}

HC.Plane.no.cota <- function(dimension, signal.values){
  signal.values = data.frame("H" = signal.values[,1], "C" = signal.values[,2], "Rho" = factor(c(rep("3.5 - 3.99", 50), "4.0")))
  
  p = ggplot(signal.values, aes(x = H, y = C, color = Rho)) + 
    geom_point(size = 2) +
    scale_shape_identity() +
    theme_few(base_size = 18, base_family = "serif")  + 
    theme(plot.title = element_text(hjust=0.5)) + 
    scale_colour_few("Dark")
  print(p)
  return(p)
}

plot.BP.analysis <- function(){
  
  plots = array(list(), 4)
  n.total = 51
  
  d = 3
  tau = 1 
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D3T3.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[1]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 3, Delay = 3")))
  
  d = 4
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D4T4.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[2]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 4, Delay = 4")))
  
  d = 5
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D5T5.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[3]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 5, Delay = 5")))
  
  d = 6
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/BP_HC_D6T6.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[4]] = HC.Plane.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 6, Delay = 6")))
  
  pdf("HCAnalysis.pdf", width = 10, height = 7) 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + 
    labs(title = "Bandt-Pompe") +
    labs(x = TeX("\\textit{Normalized Entropy}"), y = TeX("\\textit{Statistical Complexity}"))  +
    theme_igray() + theme(text = element_text(size = 16, family="Times", face="italic"), plot.title = element_text(hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  dev.off() 
}

plot.TG.analysis <- function(){
  
  plots = array(list(), 4)
  n.total = 51
  
  d = 3
  tau = 1 
  Entropy.Complexity.csv = read.csv(file="Data/TG_HC_D3T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[1]] = HC.Plane.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 3, Delay = 1")))
  
  d = 4
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/TG_HC_D4T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[2]] = HC.Plane.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 4, Delay = 1")))
  
  d = 5
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/TG_HC_D5T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[3]] = HC.Plane.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 5, Delay = 1")))
  
  d = 6
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/TG_HC_D6T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[4]] = HC.Plane.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 6, Delay = 1")))
  
  pdf("TG_T1.pdf", width = 10, height = 7) 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + 
    labs(title = "Transition Graphs") +
    labs(x = TeX("\\textit{Normalized Entropy}"), y = TeX("\\textit{Statistical Complexity}"))  +
    theme_igray() + theme(text = element_text(size = 16, family="Times", face="italic"), plot.title = element_text(hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  dev.off() 
}

plot.WPE.analysis <- function(){
  
  plots = array(list(), 4)
  n.total = 51
  
  d = 3
  tau = 1 
  Entropy.Complexity.csv = read.csv(file="Data/WPE_HC_D3T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[1]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 3, Delay = 1")))
  
  d = 4
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/WPE_HC_D4T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[2]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 4, Delay = 1")))
  
  d = 5
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/WPE_HC_D5T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[3]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 5, Delay = 1")))
  
  d = 6
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/WPE_HC_D6T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[4]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 6, Delay = 1")))
  
  pdf("WPE_T1.pdf", width = 10, height = 7) 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + 
    labs(title = "WPE") +
    labs(x = TeX("\\textit{Normalized Entropy}"), y = TeX("\\textit{Statistical Complexity}"))  +
    theme_igray() + theme(text = element_text(size = 16, family="Times", face="italic"), plot.title = element_text(hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  dev.off() 
}


plot.FGPE.analysis <- function(){
  
  plots = array(list(), 4)
  n.total = 51
  
  d = 3
  tau = 1 
  Entropy.Complexity.csv = read.csv(file="Data/FGPE_HC_D3T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[1]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 3, Delay = 1")))
  
  d = 4
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/FGPE_HC_D4T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[2]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 4, Delay = 1")))
  
  d = 5
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/FGPE_HC_D5T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[3]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 5, Delay = 1")))
  
  d = 6
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/FGPE_HC_D6T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[4]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 6, Delay = 1")))
  
  pdf("FGPE_T1.pdf", width = 10, height = 7) 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + 
    labs(title = "FGPE") +
    labs(x = TeX("\\textit{Normalized Entropy}"), y = TeX("\\textit{Statistical Complexity}"))  +
    theme_igray() + theme(text = element_text(size = 16, family="Times", face="italic"), plot.title = element_text(hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  dev.off() 
}


plot.AAPE.analysis <- function(){
  
  plots = array(list(), 4)
  n.total = 51
  
  d = 3
  tau = 1 
  Entropy.Complexity.csv = read.csv(file="Data/AAPE_HC_D3T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[1]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 3, Delay = 1")))
  
  d = 4
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/AAPE_HC_D4T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[2]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 4, Delay = 1")))
  
  d = 5
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/AAPE_HC_D5T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[3]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 5, Delay = 1")))
  
  d = 6
  tau = 1
  Entropy.Complexity.csv = read.csv(file="Data/AAPE_HC_D6T1.csv", header=TRUE, sep=",")
  Entropy.Complexity = matrix(nrow = n.total, ncol = 2)
  Entropy.Complexity[,1] = Entropy.Complexity.csv[1:n.total, 'V1']
  Entropy.Complexity[,2] = Entropy.Complexity.csv[1:n.total, 'V2']
  plots[[4]] = HC.Plane.no.cota(d, Entropy.Complexity) + ggtitle(expression(italic("Dimension = 6, Delay = 1")))
  
  pdf("AAPE_T1.pdf", width = 10, height = 7) 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            ncol = 2, nrow = 2, common.legend = TRUE, legend = "right") + 
    labs(title = "AAPE") +
    labs(x = TeX("\\textit{Normalized Entropy}"), y = TeX("\\textit{Statistical Complexity}"))  +
    theme_igray() + theme(text = element_text(size = 16, family="Times", face="italic"), plot.title = element_text(hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
  dev.off() 
}