###prepare data
###vcftools --vcf test.vcf --012 --out cs.ddRAD.gf

###cut -f2- cs.ddRAD.gf.012 | sed 's/-1/NA/g' >snp.temp
###tr -d '\t' <cs.ddRAD.gf.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
###paste <(echo "ID" | cat - cs.ddRAD.gf.012.indv) <(echo "" | cat header - snp.temp) > snp.forR
###rm header snp.temp

library(sp)
library(raster)
library(rgdal)

require(raster)
require(rgdal)


snp <- read.table("snp.forR", header = T, row.names = 1)

sample.coord <-read.table("tp.sample.location.south.txt", header=T, stringsAsFactors=F)


library(permute)
library(lattice)

library(vegan)
coord <- sample.coord[,c("longitude","latitude")]
pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
keep <- round(length(which(pcnm$value > 0))/2)
pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
write.table(pcnm.keep, "pcnm.keep", sep="\t", quote=F, row.names=F) 
pcnm.keep

clim.points <- read.table("tp.3climate", header = T, stringsAsFactors=F)

env.gf <- cbind(clim.points, pcnm.keep)

library(gradientForest)


maxLevel <- log2(0.368*nrow(env.gf)/2)
gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=1000, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
#gf<- missForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=1000, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
save.image("tp.gf.data")

###plot bar graphs depicting the importance of each spatial and climate variable.
pdf("GF_VariableImportance.pdf")
plot(gf, plot.type = "O")
dev.off()

###plot the "turnover functions" showing how allelic composition changes along the spatial or environmental gradients. 
by.importance <- names(importance(gf))

pdf("GF_TurnoverFunctions.pdf")
plot(gf, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
dev.off()

cu.sp.bio05 <- cumimp(gf, "bio_5", "Overall")
cu.sp.bio08 <- cumimp(gf, "bio_8", "Overall")
cu.sp.bio10 <- cumimp(gf, "bio_10", "Overall")
cu.sp.PCNM1 <- cumimp(gf, "pcnm.keep", "Overall")


write.table(cu.sp.bio05,file="cs.bio05.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.bio08,file="cs.bio08.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.bio10,file="cs.bio10.cumimp.PD.res",sep = " ",quote = F)
write.table(cu.sp.PCNM1,file="cs.PCNM1.cumimp.PD.res",sep = " ",quote = F)


gm.imp <- as.data.frame(t(gf$imp.rsq))

write.table(gm.imp,"tp.gf.snp.imp",sep = " ",row.names=T,col.names = T)

for (i in 1:length(gm.imp[,"bio_5"])) {
  bar <- gm.imp[i,]
  gm.imp[i,5] <- names(which.max(abs(bar[1:4]))) # gives the variable, j in gm.imp[i,j] is 3+biovars+PCNMvars+1, m and n in bar[m:n] is m=4, n= 4+vasrs(biovars + PCNMvars)
  gm.imp[i,6] <- max(abs(bar[1:4]))              # gives the correlation,j in gm.imp[i,j] is 3+biovars+PCNMvars+2, m and n in bar[m:n] is m=4, n= 4+vasrs(biovars + PCNMvars)
}

colnames(gm.imp)[5] <- "predictor"
colnames(gm.imp)[6] <- "correlation"

table(gm.imp$predictor)
gm.imp.bios <- gm.imp[which(gm.imp$predictor%in%c("bio_5","bio_8","bio_10")),]
write.table(gm.imp.bios,"gm.gf.snp.imp.4bios",sep = " ",row.names=T,col.names = T)

quants.cor <- quantile(gm.imp.bios[,"correlation"],probs = c(0.5,0.9,0.95,0.99,0.999),names = T)
quants.cor[1]
gm.imp.bios[gm.imp.bios[,"correlation"]<quants.cor[1],] <- NA

gm.imp.bios.high <- gm.imp.bios[apply(gm.imp.bios[,1:6],1,function(gm.imp.bios)!all(is.na(gm.imp.bios))),]
write.table(gm.imp.bios.high,"cs.gf.imp.snp.bios.cor.0.5.bios",quote = F,sep = " ",row.names = T,col.names = T)

