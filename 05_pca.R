#install.packages(c("poppr","vcfR",RColorBrewer,ggplot2,reshape2,cowplot,adegenet,"snowfall")iii:)
#install.packages("vcfR")
#
#install.packages(c("RColorBrewer","ggplot2","reshape2","cowplot","adegenet","snowfall"))
library("ade4")
library("adegenet")
##library("poppr")
library("vcfR")
#library(RColorBrewer)
##library("ggplot2")
##library("reshape2")
##library("cowplot")
##library("snow")
##library("snowfall")
#install.packages("snowfall")

### import vcf file###
vcf <- vcfR::read.vcfR("D:/00_caolijunamd/input.recode.vcf", verbose = FALSE) 
####converge to genlight####
genlight <- vcfR2genlight(vcf)
#pop <- as.factor(substr(genlight@ind.names,1,4))

pca2 <- prcomp(genlight, scale=T)

#write.csv(pop,file = "pop.id.csv",quote = F,row.names = F)
####define population of each sample####
#popmap <-read.table("pop.id.csv",colClasses = c("character"))# [-c(21:88,113:132,199:218),]
pop(genlight) <- pop  ### 
#pop<- as.factor(t(popmap)) ### for color assinmentw
#pop(genlight) <- as.factor(c("BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX","BJYQX","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP"))
#####################
#####simple PCA####
######################
x <- tab(genlight,NA.method="mean")
#x1 <- x[-c(21:88,113:132,199:218),]
#x1 <- x1[,which(colSums(x1)>0)]
#BJYQ <- x[21:40,]
#BJYQ2 <- prcomp(BJYQ, scale=T)
#x[2:40,2:40]
#pca2 <- prcomp(x, scale=T)
#x11(height=6, width=12, pointsize=12); par(mfrow=c(1,2)) 
mycolors <- rainbow(19) 

xlab <- paste("PC1","(",round((summary(pca2))$importance[2,1]*100,1),"%)",sep="")
ylab <- paste("PC2","(",round((summary(pca2))$importance[2,2]*100,1),"%)",sep="")

pdf(file="1pca.input.recode.vcf.pdf",width = 12,height = 8)    
plot(pca2$x[,1:2],pch=20,cex=1.8, col=mycolors[pop],xlab = xlab, ylab = ylab)
legend("topleft",ncol=4,legend=levels(pop),bty="n",col="gray32",pch=21,cex=1.5,pt.bg=mycolors)
dev.off()

pdf(file="2pca.input.recode.vcf.pdf",width = 12,height = 8)    
plot(pca2$x, type="n",xlab = xlab, ylab = ylab); text(pca2$x, rownames(pca2$x), cex=1.5, col=mycolors[pop]) 
legend("topleft",ncol=4,legend=levels(pop),bty="n",col="gray32",pch=21,cex=1.5,pt.bg=mycolors)
dev.off()
