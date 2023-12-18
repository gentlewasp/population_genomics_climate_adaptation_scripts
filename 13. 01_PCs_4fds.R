library("vegan")
library("data.table")
Neutral <- read.table("keep_field_121ind.4fds.raw",header = TRUE)
Neutral = as.data.frame(Neutral)
rownames(Neutral) = Neutral[,1] # name row names for table, the first column of table was used
Neutral = Neutral[,-(1:6)] 
Neutral <- apply(Neutral,2,function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
pca <- rda(Neutral, scale=T) # PCA in vegan uses the rda() call without any predictors

screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")
