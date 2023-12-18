#plink --vcf keep_field_121ind.recode.vcf --recode A --allow-extra-chr --out keep_field_121ind.allsnp
library("vegan")
library("data.table")

#genotype <- fread("16pops.210ind.allsnp.raw",header = TRUE)
#genotype = as.data.frame(genotype)
#rownames(genotype) = genotype[,1] # name row names for table, the first column of table was used
#genotype = genotype[,-(1:6)] 
#genotype <- apply(genotype,2,function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#
#Neutral <- read.table("16pops.210ind.raw",header = TRUE)
#Neutral = as.data.frame(Neutral)
#rownames(Neutral) = Neutral[,1] # name row names for table, the first column of table was used
#Neutral = Neutral[,-(1:6)] 
#Neutral <- apply(Neutral,2,function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#pca <- rda(Neutral, scale=T) # PCA in vegan uses the rda() call without any predictors
#
#save.image("01_fread_tp_rda_16pops.RData")
#

load("01_fread_tp_rda_16pops.RData")
PCs <- scores(pca, choices=c(1:4), display="sites", scaling=0)
PopStruct <- data.frame(Individuals = Neutral[,1], PCs)


### Table gathering all variables
Variables1 <- read.table("4bios_geo_16pops.env",header = TRUE)
## Standardization of the variables

## Table gathering all variables
Variables <- data.frame(Variables1, PopStruct[,-1])

## Full model
RDAfull <- rda(genotype ~ ., Variables)


###Variance partitioning: disentangling the drivers of genetic variation
## Full model
pRDAfull <- rda(genotype ~ PC1 + PC2 + PC3 + PC4 + long + lat +  bio05 + bio11 + bio13 + bio17,  Variables)

sink("rda.results")
RsquareAdj(pRDAfull)
anova(pRDAfull)
sink(NULL)

save.image("02_pRDAfull_tp_rda_16pops.RData")
### Pure climate model
pRDAclim <- rda(genotype ~ bio05 + bio11 + bio13 + bio17 + Condition(long + lat + PC1 + PC2 + PC3 + PC4),  Variables)

sink("rda.results", append=TRUE)
RsquareAdj(pRDAclim)
anova(pRDAclim)
sink(NULL)

save.image("03_pRDAclim_tp_rda_16pops.RData")
### Pure neutral population structure model  
pRDAstruct <- rda(genotype ~ PC1 + PC2 + PC3 + PC4 + Condition(long + lat + bio05 + bio11 + bio13 + bio17),  Variables)
sink("rda.results", append=TRUE)
RsquareAdj(pRDAstruct)
anova(pRDAstruct)
sink(NULL)

save.image("04_pRDAstruct_tp_rda_16pops.RData")

###Pure geography model 
pRDAgeog <- rda(genotype ~ long + lat + Condition(bio05 + bio11 + bio13 + bio17 + PC1 + PC2 + PC3 + PC4),  Variables)

sink("rda.results", append=TRUE)
RsquareAdj(pRDAgeog)
anova(pRDAgeog)
sink(NULL)


save.image("05_pRDAgeog_tp_rda_16pops.RData")
#load("fww_rda_all.RData")

RDA_env_full <- rda(genotype ~ bio05 + bio11 + bio13 + bio17 + PC1 + PC2 + PC3 + PC4,  Variables)

locus_scores_sample_full <- scores(RDA_env_full, choices=c(1:2), display="sites", scaling="none") # vegan references "species", here these are the loci
TAB_sample_full <- data.frame(names = row.names(locus_scores_sample_full), locus_scores_sample_full)
TAB_sample_full$pop <- c(rep("BJCY", 15), rep("BJDX", 15), rep("GDSZ", 15), rep("GDZQ", 10), rep("HNCS", 15), rep("HNS1", 15), rep("HNSY", 11), rep("LNAS", 15), rep("NMHS", 13), rep("SCCD", 15), rep("SCDY", 12), rep("SDS1", 10), rep("SDS3", 8), rep("SDSG", 13), rep("YNBN", 13), rep("YNXS", 15))
TAB_sample_full <- TAB_sample_full[order(TAB_sample_full$pop),]
TAB_var_full <- as.data.frame(scores(RDA_env_full, choices=c(1,2), display="bp")) # pull the biplot scores
## Biplot of RDA sample and variables scores
library(ggplot2)
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_sample_full, aes(x=RDA1*3, y=RDA2*3, colour = pop), size = 1.4) +
  scale_color_manual(values = rainbow(16)) +
  geom_segment(data = TAB_var_full, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_full, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var_full)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
ggsave("4.1_2_Biplot_of_RDA_sample_and_variables_scores_fullmodel.pdf",width = 8, height = 6)

save.image("06_tp_rda_16pops.RData")

####Genotype-Environment Associations: identifying loci under selection
## Conducting the genome scan using pRDA
RDA_env <- rda(genotype ~ bio05 + bio11 + bio13 + bio17 + Condition(PC1 + PC2 + PC3 + PC4),  Variables)
#choose a number of RDA axes to include when conducting the genome scan
#plot
pdf("Eigenvalues_of_constrained_axes.pdf",width=6,height=6)
screeplot(RDA_env, main="Eigenvalues of constrained axes")
dev.off()

save.image("07_GEA_tp_rda_16pops.RData")
## Formatting table for ggplot
locus_scores_sample <- scores(RDA_env, choices=c(1:2), display="sites", scaling="none") # vegan references "species", here these are the loci
TAB_sample <- data.frame(names = row.names(locus_scores_sample), locus_scores_sample)
TAB_sample$pop <- c(rep(rep("BJCY", 15), rep("BJDX", 15), rep("GDSZ", 15), rep("GDZQ", 10), rep("HNCS", 15), rep("HNS1", 15), rep("HNSY", 11), rep("LNAS", 15), rep("NMHS", 13), rep("SCCD", 15), rep("SCDY", 12), rep("SDS1", 10), rep("SDS3", 8), rep("SDSG", 13), rep("YNBN", 13), rep("YNXS", 15))
TAB_sample <- TAB_sample[order(TAB_sample$pop),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores
# Biplot of RDA sample and variables scores
library(ggplot2)
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_sample, aes(x=RDA1*3, y=RDA2*3, colour = pop), size = 1.4) +
  scale_color_manual(values = rainbow(16)) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
ggsave("4.1_2_Biplot_of_RDA_sample_and_variables_scores.pdf",width = 8, height = 6)



## Function rdadapt
source("rdadapt.R")
library(robust)
library(qvalue)
## Running the function with K = 2
rdadapt_env<-rdadapt(RDA_env, 2)

## P-values threshold after Bonferroni correction
#thres_env <- 0.1/length(rdadapt_env$p.values)
thres_env <- 0.001

## Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

## List of outlier names
outliers <- outliers[order(outliers$contig, outliers$p.value),]
write.table(outliers,"fww.rda.outlier_all",quote=F,row.names=F,col.names=F)

## Top hit outlier per contig
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])
write.table(outliers_rdadapt_env,"fww.rda.outlier_top",quote=F,row.names=F,col.names=F)

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores
## Biplot of RDA loci and variables scores
library(ggplot2)
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 0.7) +
  scale_color_manual(values = c("gray70", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
ggsave("4.1_2_Biplot_of_RDA_loci_and_variables_scores.pdf",width = 8, height = 6)
## Manhattan plot
Outliers <- rep("Neutral", length(colnames(genotype)))
Outliers[colnames(genotype)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(genotype)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(genotype)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray70", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

#ggsave("Manhattan_plot.pdf",device = cairo_pdf,width =3.15, height =2.36)
ggsave("4.1_3_Manhattan_plot.pdf",width = 16, height = 4)

### Not accounting for population structure
# Running a simple RDA model
RDA_env_unconstrained <- rda(genotype ~ bio05 + bio11 + bio13 + bio17,  Variables)
# Running the rdadapt function
rdadapt_env_unconstrained <- rdadapt(RDA_env_unconstrained, 2)
# Setting the p-value threshold 
thres_env <- 0.001
#thres_env <- 0.1/length(rdadapt_env_unconstrained$p.values)

## Identifying the outliers for the simple RDA
outliers_unconstrained <- data.frame(Loci = colnames(genotype)[which(rdadapt_env_unconstrained$p.values<thres_env)], p.value = rdadapt_env_unconstrained$p.values[which(rdadapt_env_unconstrained$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype)[which(rdadapt_env_unconstrained$p.values<thres_env)], split = ":"), function(x) x[1])))
outliers_unconstrained <- outliers_unconstrained[order(outliers_unconstrained$contig, outliers_unconstrained$p.value),]
write.table(outliers_unconstrained,"fww.rda.outlier_all_unconstrained",quote=F,row.names=F,col.names=F)

outliers_rdadapt_env_unconstrained <- as.character(outliers_unconstrained$Loci[!duplicated(outliers_unconstrained$contig)])
write.table(outliers_rdadapt_env_unconstrained,"fww.rda.outlier_top_unconstrained",quote=F,row.names=F,col.names=F)

##compared the outliers identified when accounting or not for population structure.
### For all the outliers
#library(ggVennDiagram)
list_outliers_RDA_all <- list(RDA_constrained = as.character(outliers$Loci), RDA_unconstrained = as.character(outliers_unconstrained$Loci))
#ggVennDiagram(list_outliers_RDA_all, category.names = c("partial RDA", "simple RDA"), lty="solid", color="black", size=0.2) + 
#  scale_fill_gradient2(low = "white", high = 'gray40') + guides(fill = FALSE) + theme(text = element_text(size=16, family = "Times"))

common_outliers_RDA_all <- Reduce(intersect, list_outliers_RDA_all)
write.table(common_outliers_RDA_all,"fww.rda.outlier_all_common",quote=F,row.names=F,col.names=F)

##only for the top hit locus per contig
list_outliers_RDA_top <- list(RDA_constrained = outliers_rdadapt_env, RDA_unconstrained = outliers_rdadapt_env_unconstrained)
#ggVennDiagram(list_outliers_RDA_top, category.names = c("partial RDA", "simple RDA"), lty="solid", color="black", size=0.2) + 
#  scale_fill_gradient2(low = "white", high = 'gray40') + guides(fill = FALSE) + theme(text = element_text(size=16, family = "Times"))

common_outliers_RDA_top <- Reduce(intersect, list_outliers_RDA_top)
write.table(common_outliers_RDA_top,"fww.rda.outlier_top_common",quote=F,row.names=F,col.names=F)

save.image("08_outliers_tp_rda_16pops.RData")

### Adaptively enriched RDA
RDA_outliers <- rda(genotype[,common_outliers_RDA_top] ~ bio05 + bio11 + bio13 + bio17,  Variables)

## RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
ggsave("5.1_1_Adaptively_enriched_genetic_space_plot.png",width = 10, height = 6, dpi = 600)

save.image("09_outliers_tp_rda_16pops.RData")

