library(vioplot)
xi <- read.table("trace_tau.out",header = TRUE)

plot(1,xlim = c(0.75,6.25),ylim = c(0,1),
  ylab = expression(Branch-specific~sex-ratio~(xi)),
  xlab = "Branch",type = "n",cex.lab = 1.25)
vioplot(xi$deme_1,xi$deme_2,xi$deme_3,xi$deme_4,col = "lightblue",add = TRUE)

mean(xi$deme_1)
mean(xi$deme_2)
mean(xi$deme_3)
mean(xi$deme_4)

source("KimTree.R")
draw.tree(tree_file = "../cs.tre2",n = 3,
summary_tau_file = "summary_tau.out",
leafs = TRUE,nodes = FALSE,edges = FALSE)

