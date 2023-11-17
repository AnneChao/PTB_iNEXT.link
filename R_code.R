## Install/import packages from CRAN
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(phytools)
library(bipartite)
library(patchwork)


## Install the following packages from Anne Chao's github
library(devtools)

install_github('AnneChao/iNEXT.3D')
install_github('AnneChao/iNEXT.4steps')
install_github('AnneChao/iNEXT.beta3D')
install_github('AnneChao/iNEXT.link')

## Import the package 'iNEXT.link'
library(iNEXT.link)   ## Only functions in the package 'iNEXT.link' are used


## Because this code runs a random bootstrapping process in the background, with default = 100 replications/loops, 
## to estimate standard error, the output (for s.e. and 95% confidence intervals) will vary very slightly each time 
## you enter the same data. Also, the bootstrap procedure is VERY time-consuming; please change the argument 
## nboot = 100 (default) in any function to nboot = 0 (to skip the bootstrap process), or nboot = 20 
## (to roughly browse the resulting confidence intervals).  


## ================= Figure 1: A simple hypothetical bipartite network ================== ##
data = list("Demo" = matrix(c(6, 0, 2, 
                              5, 3, 0, 
                              1, 4, 0, 
                              1, 0, 8), ncol = 4))


rownames(data$Demo) = c("a", "b", "c")
colnames(data$Demo) = c("1", "2", "3", "4")


sequence <- list(seq.higher = colnames(data$Demo), 
                 seq.lower  = rownames(data$Demo))


png("Figure 1.png", width = 3600, height = 3000, res = 1000)
plotweb(sortweb(data$Demo, sort.order = "seq", sequence = sequence), 
        labsize = 2, col.low = 'blue', col.high = 'red', method = "normal",
        text.low.col = "blue", text.high.col = "red")
dev.off()



## ================= Figure 2: Taxonomic network diversity profiles  =========== ##
qD <- function(x, q){
  
  p <- x[x > 0] / sum(x)
  
  Sub <- function(q) {
    if (q == 0) sum(p > 0)
    
    else if (q == 1) exp( -sum(p * log(p)) )
    
    else exp(1 / (1 - q) * log( sum(p^q) ) )
  }
  
  sapply(q, Sub)
}

site1 = rep(30/8, 8)
site2 = c(1, 2, 3, 4, 4, 5, 5, 6)
site3 = c(1, 1, 2, 3, 4, 5, 6, 8)
site4 = c(1, 1, 1, 1, 1, 1, 1, 23)

q = seq(0, 3, 0.1)

pdf("Figure 2.pdf", width = 6, height = 5)

plot(q, qD(site1, q), type = 'l', lty = 1, col = 'black',
     ylim = c(1, 8.6), yaxt = "n", 
     xlab = substitute(paste(bold('Order '), bolditalic(q), sep = '')), 
     ylab = 'Taxonomic network diversity',
     font.lab = 2, lwd = 2)

axis(side = 2, at = seq(1, 9, 1), labels = seq(1, 9, 1))
lines(q, qD(site2, q), col = 'red',          lty = 'longdash', lwd = 2)
lines(q, qD(site3, q), col = 'springgreen4', lty = '1342',     lwd = 2)
lines(q, qD(site4, q), col = 'blue',         lty = 3,          lwd = 2)

text(2.6, 8.3, substitute(paste( bold('Completely even'))),   col = 'black',        font.lab = 2)
text(2.6, 7,   substitute(paste( bold('Slightly uneven'))),   col = 'red',          font.lab = 2)
text(2.6, 5,   substitute(paste( bold('Moderately uneven'))), col = 'springgreen4', font.lab = 2)
text(2.7, 2,   substitute(paste( bold('Highly uneven'))),     col = 'blue',         font.lab = 2)

dev.off()






## ============================ load data for Figures 5 ~ 9 ============================ ##
complete_data = read.csv("Data_tree-beetle_interaction_frequency.csv")
beetles_col_tree = read.tree("Data_phylo_tree.txt")
beetles_col_distM = read.table("Data_distance_matrix.txt")


## transform data into the format of 'iNEXT.link'
woNet_data =  list(G = 1, N = 1, O = 1)
for(t in 1:3) {
  
  ## select plot A only
  tem_data = filter(complete_data, treatment == c("G","N","O")[t] & plot == 'A')
  
  ## select pool data
  # tem_data = filter(complete_data,treatment == c("G","N","O")[t])
  
  tem_data[is.na(tem_data)] = 0
  tem_data = tem_data[,-1]
  tem_data1 = acast(tem_data, tree_species ~ species, sum, value.var = 'abundance')
  
  woNet_data[[t]] = tem_data1
  
}

names(woNet_data) = c('Closed', 'Net', 'Open' )
woNet_data = woNet_data[c(1,3)]    ## selece site 'Closed' and 'Open'


## =============== Figure 5: Plots of sample network completeness curves ============== ##
output.SC = Completeness.link(woNet_data, nboot = 100)

ggsave(filename = 'Figure 5.pdf', 
       
       ggCompleteness.link(output.SC) + 
         scale_colour_manual(values = c('blue', 'red')) +
         scale_fill_manual(values = c('blue', 'red')),
       
       width = 5, height = 4, dpi = 1000)



## ===================================== Figure 6 & 7 ===================================== ##

## ========= Figures 6(a1) & 7(a): Computing taxonomic network diversity using functions from iNEXT.link ================== ##
output.iNEXT.TD = iNEXT.link(data = woNet_data, diversity = 'TD', q = c(0,1,2), nboot = 100)
output.est.TD = estimateD.link(woNet_data, q = c(1,2), base = "coverage", level = c(0.95,1), nboot = 100)

output.iNEXT.TD$TDiNextEst$coverage_based = rbind(output.iNEXT.TD$TDiNextEst$coverage_based, 
                                                  output.est.TD %>% select(-"s.e."))

ggiNEXT.TD = ggiNEXT.link(output.iNEXT.TD, facet.var = "Order.q")

size.TD = ggiNEXT.TD[[1]] + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ggtitle("(a1) Taxonomic network diversity") + ylab("Taxonomic network diversity") +
  theme(legend.position = "none", plot.title = element_text(color = "blue", size = 25))

cov.TD = ggiNEXT.TD[[3]] + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ggtitle("(a) Taxonomic network diversity") + ylab("Taxonomic network diversity") + 
  theme(legend.position = "none", plot.title = element_text(color = "blue", size = 25)) + 
  coord_cartesian(xlim = c(0.5,1))


## =========== Figures 6(b1) & 7(b): Computing phylogenetic network diversity using functions from iNEXT.link ======================= ##
output.iNEXT.PD = iNEXT.link(data = woNet_data, diversity = 'PD', q = c(0,1,2), nboot = 100, col.tree = beetles_col_tree)
output.est.PD = estimateD.link(woNet_data, diversity = "PD", q = c(1,2), base = "coverage", level = c(0.95,1), nboot = 100, col.tree = beetles_col_tree)

output.iNEXT.PD$PDiNextEst$coverage_based = rbind(output.iNEXT.PD$PDiNextEst$coverage_based, 
                                                  output.est.PD %>% select(-"s.e."))

ggiNEXT.PD = ggiNEXT.link(output.iNEXT.PD, facet.var = "Order.q")

size.PD = ggiNEXT.PD[[1]] + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ggtitle("(b1) Mean phylogenetic network diversity") + ylab("Mean phylogenetic network diversity") + 
  theme(legend.position = "none", plot.title = element_text(color = "blue", size = 25)) + 
  facet_wrap(~ paste("q =", Order.q), nrow = 1)

cov.PD = ggiNEXT.PD[[3]] + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ggtitle("(b) Mean phylogenetic network diversity") + ylab("Mean phylogenetic network diversity") +
  theme(legend.position = "none", plot.title = element_text(color = "blue", size = 25)) + 
  facet_wrap(~ paste("q =", Order.q), nrow = 1) + 
  coord_cartesian(xlim = c(0.5,1))


## =============== Figures 6(c1) & 7(c): Computing functional network diversity using functions from iNEXT.link ================ ##
output.iNEXT.FD = iNEXT.link(data = woNet_data, diversity = 'FD', q = c(0,1,2), nboot = 100, col.distM = beetles_col_distM)
output.est.FD = estimateD.link(woNet_data, diversity = "FD", q = c(1,2), base = "coverage", level = c(0.95,1), nboot = 100, col.distM = beetles_col_distM)

output.iNEXT.FD$FDiNextEst$coverage_based = rbind(output.iNEXT.FD$FDiNextEst$coverage_based, 
                                                  output.est.FD %>% select(-'s.e.'))

ggiNEXT.FD = ggiNEXT.link(output.iNEXT.FD, facet.var = "Order.q")

size.FD = ggiNEXT.FD[[1]] + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ggtitle("(c1) Functional network diversity") + ylab("Functional network diversity") +
  theme(plot.title = element_text(color = "blue", size = 25))

cov.FD = ggiNEXT.FD[[3]] + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ggtitle("(c) Functional network diversity") + ylab("Functional network diversity") +
  theme(plot.title = element_text(color = "blue", size = 25)) + 
  coord_cartesian(xlim = c(0.5,1))



Figure_7 = cov.TD / cov.PD / cov.FD

ggsave(filename = 'Figure 7.pdf', Figure_7, width = 10, height = 16, dpi = 1000)          ## Figure 7 ##


## ================ Figure 6(a2): compute asymptotic & observed taxonomic network diversity =========================== ##
AO.TD = ObsAsy.link(data = woNet_data, diversity = 'TD', q = seq(0, 2, 0.2), nboot = 100)

ggAO.TD = ggObsAsy.link(AO.TD) + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ylab("Taxonomic network diversity") +
  ggtitle("(a2) Asy. and obs. profiles") + 
  theme(legend.position = "none",
        plot.title = element_text(color = "blue", size = 25),
        text = element_text(size = 16))


## ================== Figure 6(b2): compute asymptotic & observed phylogenetic network diversity =========================== ##
AO.PD = ObsAsy.link(data = woNet_data, diversity = 'PD', q = seq(0, 2, 0.2), nboot = 100, col.tree = beetles_col_tree)

ggAO.PD = ggObsAsy.link(AO.PD) + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ylab("Mean phylogenetic network diversity") +
  ggtitle("(b2) Asy. and obs. profiles") + 
  theme(legend.position = "none",
        plot.title = element_text(color = "blue", size = 25),
        text = element_text(size = 16))


## =============== Figure 6(c2): compute asymptotic & observed functional network diversity ===================== ##
AO.FD = ObsAsy.link(data = woNet_data, diversity = 'FD', q = seq(0, 2, 0.2), nboot = 100, col.distM = beetles_col_distM)

ggAO.FD = ggObsAsy.link(AO.FD) + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ylab("Functional network diversity") +
  ggtitle("(c2) Asy. and obs. profiles") + 
  theme(plot.title = element_text(color = "blue", size = 25),
        text = element_text(size = 16))


Figure_6.TD = (size.TD + coord_cartesian(ylim = c(0,max(AO.TD$qTD.UCL))) ) + (ggAO.TD + coord_cartesian(ylim = c(0,max(AO.TD$qTD.UCL))) ) + plot_layout(widths = c(2, 1))
Figure_6.PD = (size.PD + coord_cartesian(ylim = c(0,max(AO.PD$qPD.UCL))) ) + (ggAO.PD + coord_cartesian(ylim = c(0,max(AO.PD$qPD.UCL))) ) + plot_layout(widths = c(2, 1))
Figure_6.FD = (size.FD + coord_cartesian(ylim = c(0,64)) ) + (ggAO.FD + coord_cartesian(ylim = c(0,64)) ) + plot_layout(widths = c(2, 1))
# Figure_6.FD = (size.FD + coord_cartesian(ylim = c(0,100)) ) + (ggAO.FD + coord_cartesian(ylim = c(0,100)) ) + plot_layout(widths = c(2, 1))


ggsave(filename = 'Figure 6.pdf', Figure_6.TD / Figure_6.PD / Figure_6.FD, 
       width = 16, height = 16, dpi = 1000)                                   ## Figure 6 ##



## ==================== Figure 8: compute two specialization measures ================== ##
output.Spec = Spec.link(woNet_data, E.class = c(1,3), nboot = 100)

Figure_8 = ggSpec.link(output.Spec) + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  theme(strip.text.x = element_text(size = 12, colour = "black", face = "bold"))

ggsave(filename = 'Figure 8.pdf', Figure_8, 
       width = 6, height = 4, dpi = 1000)


## ===================================== Figure in Appendix S4 ===================================== ##

S4.function = function(coverage) {
  beetles_data = list(Ap = list(G=1,N=1,O=1), Bp = list(G=1,N=1,O=1),
                      Cp = list(G=1,N=1,O=1), Dp = list(G=1,N=1,O=1),
                      Ep = list(G=1,N=1,O=1), Fp = list(G=1,N=1,O=1))
  
  pp = c("A","B","C","D","E","F")
  tt = c("G","N","O")
  
  for(p in 1:6){
    for(t in 1:3){
      tem = filter(complete_data,treatment == tt[t],plot == pp[p])
      tem_data = tem[,c(9,7)]
      tem_data = tem_data %>% group_by(inte) %>% summarise_all(sum) %>% as.data.frame() 
      rownames_tem = tem_data[,1]
      tem_data = tem_data[,-1]
      tem_data = as.matrix(tem_data)
      rownames(tem_data) = rownames_tem
      colnames(tem_data) = "."
      ttem = iNEXT.link:::long_to_wide(tem_data)
      
      beetles_data[[p]][[t]] = ttem
    }
  }
  
  out = data.frame()
  
  for(i in 1:6){
    names(beetles_data[[i]]) = c('Closed','Net','Open')
    woNet_data_tem = list(Closed = beetles_data[[i]][[1]], Open = beetles_data[[i]][[3]])
    
    output.Spec = Spec.link(woNet_data_tem, q = c(1,2), E.class = 1, C = coverage, nboot = 0)[[1]]
    output.Spec$Site = paste("plot", c("A", "B", "C", "D", "E", "F")[i])
    output.Spec$Order.q = paste("q =", output.Spec$Order.q, ",  C = ", coverage*100, "%")
    
    H = c(H2fun(woNet_data_tem$Closed)[1],
          H2fun(woNet_data_tem$Open)[1])
    
    output.H2 = data.frame(Order.q = "H2'", 
                           Specialization = H, 
                           Site = paste("plot", c("A", "B", "C", "D", "E", "F")[i]), 
                           Network = c("Closed","Open"))
    
    out = rbind(out, output.Spec[, c('Order.q', 'Specialization', 'Site', 'Network')], 
                output.H2)
  }
  
  return(out)
}



## Figure in Appendix S4, coverage = 90%
output.S4.0.9 = S4.function(0.9)

Figure_S4_0.9 = ggplot(output.S4.0.9, aes(x = Network, y = Specialization, col = Site, group = Site)) +
  geom_line(size = 0.4, lty = 1) +
  geom_point(aes(shape = Site), size = 2) +
  labs(x = "Habitat") +
  theme_bw() +
  facet_grid(. ~ Order.q) +
  scale_colour_manual(values = c('#F8766D', '#BB9D00', 'purple',
                                 '#00C0B8', 'gray60', '#E76BF3'))


## Figure in Appendix S4, coverage = 95%
output.S4.0.95 = S4.function(0.95)

Figure_S4_0.95 = ggplot(output.S4.0.95, aes(x = Network, y = Specialization, col = Site, group = Site)) +
  geom_line(size = 0.4, lty = 1) +
  geom_point(aes(shape = Site), size = 2) +
  labs(x = "Habitat") +
  theme_bw() +
  facet_grid(. ~ Order.q) +
  scale_colour_manual(values = c('#F8766D', '#BB9D00', 'purple',
                                 '#00C0B8', 'gray60', '#E76BF3'))


ggsave(filename = 'Figure S4.pdf', Figure_S4_0.9 / Figure_S4_0.95, 
       width = 6, height = 6, dpi = 1000)                           


