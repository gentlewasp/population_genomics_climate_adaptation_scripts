library(ggplot2)
data <- read.table("19pops_DP3.GQ20.allele.imiss",head=F)
colnames(data)[1] <- "id"
colnames(data)[2] <- "missing"
data[,3] <- data[,2]/12332601
colnames(data)[3] <- "rate"
str(data)
dp.thresh=quantile(data[,3],probs=0.85,na.rm=T)
dp.thresh
dp.thresh=quantile(data[,3],probs=0.9,na.rm=T)
dp.thresh

data[data[,3]>=dp.thresh,]
p <- ggplot(data,aes(reorder(id,-rate),rate))+geom_bar(stat = 'identity')

pdf("19pops_sample_missing_1.pdf",width = 30, height = 6)
p
dev.off()
