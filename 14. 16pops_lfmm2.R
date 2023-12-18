library(LEA)
library(lfmm)
library(qqman)
library(qvalue)
library(dplyr)

LEA::vcf2geno("16pops_210samples.vcf") #convert vcf to geno
geno <- read.geno("16pops_210samples.geno")
env <- read.table("4factors_16pops.bios",head=T)

vcf2lfmm("16pops_210samples.vcf", force=T) ###convert vcf to lfmm
pca = pca("16pops_210samples.lfmm", scale=T) 
tw = tracy.widom(pca) # Perfom Tracy-Widom tests on all eigenvalues.
tw$pvalues[1:10]

pdf(file="tp.pca.explained.by.each.component_16pops_210samples.pdf", width=8, height=8)
plot(tw$percentage,xlim = c(1,15))  # plot the percentage of variance explained by each component
dev.off()

mod.lfmm <- lfmm_ridge(geno,env, K = 5)
pv <- lfmm_test(geno,env,lfmm = mod.lfmm,calibrate = "gif")
pvalues <- pv$calibrated.pvalue 

pdf("tp_qqplot_16pops_210samples.pdf",width=10,height=10)
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
dev.off()


q.bios <- matrix(nrow =nrow(pvalues), ncol =ncol(pvalues)) ####number of SNP, number of bios
for (i in 1:6){
  q.bios[,i] <- qvalue(pvalues[,i])$qvalues
}

q.bios[q.bios==0] <- min(q.bios[q.bios>0]) 
#pvalues[is.na(pvalues)%>%which] <- 0.9999 #填充missing
#pvalues[pvalues==0] <- min(pvalues[pvalues>0]) #替换0值


data <- cbind(read.table("16pops_210samples.vcfsnp")[,c(1,2)],q.bios)
data$SNP <- seq(1, nrow(data), 1)
colnames(data)[1] <- "CHROM"
colnames(data)[2] <- "POS"
colnames(data)[c(3:(ncol(q.bios)+2))] <- colnames(env)
threshold=0.01

for (i in 1:6){
data1 <- data[,c(1,2,i+2)]
outliers <- data1[data1[,3]<=threshold,]
write.table(outliers,paste("tp_lfmm_outliers_0.01_",colnames(data1[3]),sep=""),sep = "\t",
row.names = F,col.names=F,quote = F)
          }

##pdf("tp_manhattan_104samples_1-7.pdf",width=16,height=8)
jpeg("tp_manhattan_16pops_210samples_1-6.jpg",width = 3000,height = 3000)
par(mfrow = c(7,1))
for (i in 1:6){
data1 <- data[,c(1,2,i+2,9)]
maxQ <- max(-log(data1[,3],10),na.rm=T)
maxQ
outliers <- data1[data1[,3]<=threshold,]
outliers_list <- outliers$SNP
manhattan(data1,chr = "CHROM",bp = "POS",p = colnames(data1[3]),snp = "SNP",logp=T,highlight = outliers_list,
          suggestiveline = F,ylab=colnames(data1[3]),ylim=c(0,maxQ))
          }
dev.off()


#jpeg("tp_manhattan_16pops_210samples_8-14.jpg",width = 3000,height = 3000)
#par(mfrow = c(7,1))
#for (i in 8:14){
#data1 <- data[,c(1,2,i+2,24)]
#maxQ <- max(-log(data1[,3],10),na.rm=T)
#maxQ
#outliers <- data1[data1[,3]<=threshold,]
#outliers_list <- outliers$SNP
#manhattan(data1,chr = "CHROM",bp = "POS",p = colnames(data1[3]),snp = "SNP",logp=T,highlight = outliers_list,
#          suggestiveline = F,ylab=colnames(data1[3]),ylim=c(0,maxQ))
#          }
#dev.off()
#
#jpeg("tp_manhattan_16pops_210samples_15-21.jpg",width = 3000,height = 3000)
#par(mfrow = c(7,1))
#for (i in 15:21){
#data1 <- data[,c(1,2,i+2,24)]
#maxQ <- max(-log(data1[,3],10),na.rm=T)
#maxQ
#outliers <- data1[data1[,3]<=threshold,]
#outliers_list <- outliers$SNP
#manhattan(data1,chr = "CHROM",bp = "POS",p = colnames(data1[3]),snp = "SNP",logp=T,highlight = outliers_list,
#          suggestiveline = F,ylab=colnames(data1[3]),ylim=c(0,maxQ))
#          }
#dev.off()



#manhattan one by one
threshold=0.01
jpeg("tp_manhattan_16pops_210samples_CTmax_mean_qvalue.jpg",width = 2000,height = 300)
data1 <- data[,c(1,2,3,9)]
maxQ <- max(-log(data1[,3],10),na.rm=T)
maxQ
outliers <- data1[data1[,3]<=threshold,]
outliers_list <- outliers$SNP
manhattan(data1,chr = "CHROM",bp = "POS",p = colnames(data1[3]),snp = "SNP",logp=T,highlight = outliers_list,
          suggestiveline = F,ylab=colnames(data1[3]),ylim=c(0,maxQ))
dev.off()
#pdf("tp_manhattan_104samples_CTmax.pdf",width=16,height=3)
#jpeg("tp_manhattan_104samples_CTmax.jpg",width = 2000,height = 300)
#data1 <- data[,c(1,2,3,8)]
#maxp <- max(-log(data1[,3],10),na.rm=T)
#maxp
#outliers <- data1[data1[,3]<=threshold,]
#outliers_list <- outliers$SNP
#manhattan(data1,chr = "CHROM",bp = "POS",p = "pvalues",snp = "SNP",logp=T,highlight = outliers_list,
#          suggestiveline = F,ylab="-LogP",ylim=c(0,20))
#dev.off()
#
