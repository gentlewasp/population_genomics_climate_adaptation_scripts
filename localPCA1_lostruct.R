#bcftools view ExAC.vcf -Oz -o ExAC.vcf.gz
#bcftools index ExAC.vcf.gz
library(vcfR)
library(adegenet)
library(lostruct)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)

window_size <- 100
k_kept <- 40
max_distance_between_outliers <- 100
vcf.file <- "/public1/home/caolijun/caolijun/BGI20200423/GM/SNP/vcffiles/gm_31sample_all_mac1_hwee7_chr9.recode.vcf.gz"
sites <- vcf_positions(vcf.file)
win.fn.snp <- vcf_windower(vcf.file, size=window_size, type="snp", sites=sites)
system.time( snp.pca <- eigen_windows(win.fn.snp,k=2, mc.cores=4) )
system.time( pcdist <- pc_dist( snp.pca ) )
pcdist_na <- which(is.na(pcdist), TRUE)
 na.inds <- is.na(pcdist[,1])
  if (sum(na.inds) == length(na.inds)){
    na.inds <- is.na(pcdist[,2])
  }
mds <- cmdscale( pcdist[!na.inds,!na.inds], eig=TRUE, k=k_kept )
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.regions <- region(win.fn.snp)()
win.regions$n <- 1:nrow(win.regions)
win.regions <- win.regions[!na.inds,]
win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions
for (k in 1:k_kept){
    str_pad(k, 2, pad = "0")
    name = paste("mds",str_pad(k, 2, pad = "0"),sep="")
    win.regions$tmp <- "NA"
    win.regions <- win.regions %>% rename(!!name := tmp)
  }
  #Add the MDS coordinates to each window.
  for (i in 1:k_kept){
    j = i + 5
    win.regions[,j] <- mds.coords[,i]
  }

pdf("gm_chr9.win100.MDSplots_2.pdf",height=10,width=25)
  print(
    win.regions %>%
      gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-38)]) %>%
      ggplot(.,aes(x=mid,y=value)) + geom_point() + facet_grid(mds~.,scales = "free") +
      theme_bw()
  )
  dev.off()
  saveRDS(win.regions, file = "gm_chr9.win100.windows.rds")
