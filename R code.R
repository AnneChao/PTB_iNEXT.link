## library packages from CRAN
library(patchwork)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(bipartite)
library(ape)


## install packages from Anne Chao's github
library(devtools)

install_github('AnneChao/iNEXT.3D')
install_github('AnneChao/iNEXT.4steps')
install_github('AnneChao/iNEXT.Beta3D')
install_github('AnneChao/iNEXT.link')

library(iNEXT.link)   ## Here we only use functions in the package 'iNEXT.link'


## load data
complete_data = read.csv("Data tree-beetle interaction frequency.csv")
beetles_col_tree = read.tree("Data phylo_tree.txt")
beetles_col_distM = read.table("Data distance.txt")


## transform data into the format of 'iNEXT.link'
woNet_data =  list(G = 1, N = 1, O = 1)
for(t in 1:3) {
  
  ## select pool data or plot A only
  tem_data = filter(complete_data, treatment == c("G","N","O")[t] & plot == 'A')
  # tem_data = filter(complete_data,treatment == c("G","N","O")[t])
  
  tem_data[is.na(tem_data)] = 0
  tem_data = tem_data[,-1]
  tem_data1 = acast(tem_data, tree_species ~ species, sum, value.var = 'abundance')
  
  woNet_data[[t]] = tem_data1
  
}

names(woNet_data) = c('Closed', 'Net', 'Open' )
woNet_data = woNet_data[c(1,3)]    ## selece site 'Closed' and 'Open'


## ===================================== Figure 1 ===================================== ##
output.SC = Completeness.link(woNet_data, nboot = 100)

ggsave(filename = 'Figure 1.png', 
       
       ggCompleteness.link(output.SC) + 
         scale_colour_manual(values = c('blue', 'red')) +
         scale_fill_manual(values = c('blue', 'red')),
       
       width = 5, height = 4, dpi = 300)



## ===================================== Figure 2 & 3 ===================================== ##

## =========================== iNEXT of Taxonomic diversity =========================== ##
output.iNEXT.TD = iNEXT.link(data = woNet_data, diversity = 'TD', q = c(0,1,2), nboot = 100)
output.est.TD = estimateD.link(woNet_data, q = c(1,2), base = "coverage", level = c(0.95,1), nboot = 100)

output.iNEXT.TD$iNextEst$coverage_based = rbind(output.iNEXT.TD$iNextEst$coverage_based, 
                                                output.est.TD %>% select(-"s.e."))

ggiNEXT.TD = ggiNEXT.link(output.iNEXT.TD, diversity = 'TD', facet.var = "Order.q")

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


## =========================== iNEXT of Phylogenetic diversity =========================== ##
output.iNEXT.PD = iNEXT.link(data = woNet_data, diversity = 'PD', q = c(0,1,2), nboot = 100, col.tree = beetles_col_tree)
output.est.PD = estimateD.link(woNet_data, diversity = "PD", q = c(1,2), base = "coverage", level = c(0.95,1), nboot = 100, col.tree = beetles_col_tree)

output.iNEXT.PD$PDiNextEst$coverage_based = rbind(output.iNEXT.PD$PDiNextEst$coverage_based, 
                                                  output.est.PD %>% select(-"s.e."))

ggiNEXT.PD = ggiNEXT.link(output.iNEXT.PD, diversity = 'PD', facet.var = "Order.q")

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


## =========================== iNEXT of Functional diversity =========================== ##
output.iNEXT.FD = iNEXT.link(data = woNet_data, diversity = 'FD', q = c(0,1,2), nboot = 100, col.distM = beetles_col_distM)
output.est.FD = estimateD.link(woNet_data, diversity = "FD", q = c(1,2), base = "coverage", level = c(0.95,1), nboot = 100, col.distM = beetles_col_distM)

output.iNEXT.FD$AUCiNextEst$coverage_based = rbind(output.iNEXT.FD$AUCiNextEst$coverage_based, 
                                                   output.est.FD %>% select(-'s.e.'))

ggiNEXT.FD = ggiNEXT.link(output.iNEXT.FD, diversity = 'FD', facet.var = "Order.q")

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



Figure_3 = cov.TD / cov.PD / cov.FD

ggsave(filename = 'Figure 3.png', Figure_3, width = 10, height = 16, dpi = 300)          ## Figure 3 ##


## =========================== Asymptotic & Observed Taxonomic diversity =========================== ##
AO.TD = AO.link(data = woNet_data, diversity = 'TD', q = seq(0, 2, 0.2), nboot = 100)

ggAO.TD = ggAO.link(AO.TD, diversity = 'TD') + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ylab("Taxonomic network diversity") +
  ggtitle("(a2) Asy. and obs. profiles") + 
  theme(legend.position = "none",
        plot.title = element_text(color = "blue", size = 25),
        text = element_text(size = 16))


## =========================== Asymptotic & Observed Phylogenetic diversity =========================== ##
AO.PD = AO.link(data = woNet_data, diversity = 'PD', q = seq(0, 2, 0.2), nboot = 100, col.tree = beetles_col_tree)

ggAO.PD = ggAO.link(AO.PD, diversity = 'PD') + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ylab("Mean phylogenetic network diversity") +
  ggtitle("(b2) Asy. and obs. profiles") + 
  theme(legend.position = "none",
        plot.title = element_text(color = "blue", size = 25),
        text = element_text(size = 16))


## =========================== Asymptotic & Observed Functional diversity =========================== ##
AO.FD = AO.link(data = woNet_data, diversity = 'FD', q = seq(0, 2, 0.2), nboot = 100, col.distM = beetles_col_distM)

ggAO.FD = ggAO.link(AO.FD, diversity = 'FD') + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  ylab("Functional network diversity") +
  ggtitle("(c2) Asy. and obs. profiles") + 
  theme(plot.title = element_text(color = "blue", size = 25),
        text = element_text(size = 16))


Figure_2.TD = (size.TD + coord_cartesian(ylim = c(0,max(AO.TD$qD.UCL))) ) + (ggAO.TD + coord_cartesian(ylim = c(0,max(AO.TD$qD.UCL))) ) + plot_layout(widths = c(2, 1))
Figure_2.PD = (size.PD + coord_cartesian(ylim = c(0,max(AO.PD$qD.UCL))) ) + (ggAO.PD + coord_cartesian(ylim = c(0,max(AO.PD$qD.UCL))) ) + plot_layout(widths = c(2, 1))
Figure_2.FD = (size.FD + coord_cartesian(ylim = c(0,64)) ) + (ggAO.FD + coord_cartesian(ylim = c(0,64)) ) + plot_layout(widths = c(2, 1))
# Figure_2.FD = (size.FD + coord_cartesian(ylim = c(0,100)) ) + (ggAO.FD + coord_cartesian(ylim = c(0,100)) ) + plot_layout(widths = c(2, 1))


ggsave(filename = 'Figure 2.png', Figure_2.TD / Figure_2.PD / Figure_2.FD, 
       width = 16, height = 16, dpi = 300)                                   ## Figure 2 ##



## ===================================== Figure 4 ===================================== ##
output.Spec = Spec.link(woNet_data, E.class = c(1,3), nboot = 100)
Figure_4 = ggSpec.link(output.Spec) + 
  scale_colour_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) + 
  theme(strip.text.x = element_text(size = 12, colour = "black", face = "bold"))

ggsave(filename = 'Figure 4.png', Figure_4, 
       width = 6, height = 4, dpi = 300)


## ===================================== Figure 5 ===================================== ##

Fig5.function = function(coverage) {
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
    
    output.Spec = Spec.link(woNet_data_tem, q = c(1,2), E.class = 1, C = coverage, nboot = 0)[[2]]
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



## Figure 5. coverage = 90%
output.Fig5_0.9 = Fig5.function(0.9)

Figure_5_0.9 = ggplot(output.Fig5_0.9, aes(x = Network, y = Specialization, col = Site, group = Site)) +
  geom_line(size = 0.4, lty = 1) +
  geom_point(aes(shape = Site), size = 2) +
  labs(x = "Habitat") +
  theme_bw() +
  facet_grid(. ~ Order.q) +
  scale_colour_manual(values = c('#F8766D', '#BB9D00', 'purple',
                                 '#00C0B8', 'gray60', '#E76BF3'))


## Figure 5. coverage = 95%
output.Fig5.0.95 = Fig5.function(0.95)

Figure_5_0.95 = ggplot(output.Fig5.0.95, aes(x = Network, y = Specialization, col = Site, group = Site)) +
  geom_line(size = 0.4, lty = 1) +
  geom_point(aes(shape = Site), size = 2) +
  labs(x = "Habitat") +
  theme_bw() +
  facet_grid(. ~ Order.q) +
  scale_colour_manual(values = c('#F8766D', '#BB9D00', 'purple',
                                 '#00C0B8', 'gray60', '#E76BF3'))


ggsave(filename = 'Figure 5.png', Figure_5_0.9 / Figure_5_0.95, 
       width = 6, height = 6, dpi = 300)                            ## Figure 5


