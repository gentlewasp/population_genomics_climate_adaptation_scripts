library(purrr)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
################################################
# 导入自己构建的 OrgDb
################################################
#install.packages("org.Tpalmi.eg.db", repos=NULL, type="sources")
library(org.Tpalmi.eg.db)
columns(org.Tpalmi.eg.db)
keytypes(org.Tpalmi.eg.db)
keys(org.Tpalmi.eg.db, keytype="Pathway") %>% head()

# 导入需要进行富集分析的基因列表，并转换为向量
#########################################################################################
DD<-"1_2_3_BJDX_vs_SJY4.gene"
DEGs<- read.table(DD, header=F, sep = "\t")
gene_list <- DEGs[,1]

################################################
# 从 OrgDB 提取 Pathway 和基因的对应关系
################################################
pathway2gene <- AnnotationDbi::select(org.Tpalmi.eg.db, 
                                       keys = keys(org.Tpalmi.eg.db), 
                                       columns = c("Pathway","Ko")) %>%
   na.omit() %>%
   dplyr::select(Pathway, GID)
   
################################################
# 导入 Pathway 与名称对应关系
################################################
load("kegg_info.RData")

#KEGG pathway 富集
ekp <- enricher(gene_list, 
                TERM2GENE = pathway2gene, 
                TERM2NAME = pathway2name, 
                pvalueCutoff = 0.01, 
                qvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                minGSSize = 1)
                
ekp_results <- as.data.frame(ekp)
write.table(ekp_results, file = "F0_F4_ekp.txt", quote = F)                    ###让保存的字符串不用“”引起来

pdf(file = "F0_F4_ekp.pdf")                   
barplot(ekp, showCategory=20,color="pvalue",
        font.size=10)
dotplot(ekp)
x2 <- pairwise_termsim(ekp)
emapplot(x2)
cnetplot(ekp, showCategory = 5)
dev.off()                   

#########################################################################################
#GO 分析
#########################################################################################

ego <- enrichGO(gene = gene_list,                       #差异基因 vector 
                keyType = "GID",                                   #差异基因的 ID 类型，需要是 OrgDb 支持的 
                OrgDb = org.Tpalmi.eg.db,                               #对应的OrgDb 
                ont = "BP",                                             #GO 分类名称，CC BP MF 
                pvalueCutoff = 0.01,     #Pvalue 阈值 （pvalue=1指输出所有结果，pvalue=0.05指输出符合要求的结果） 
                qvalueCutoff = 0.05,     #qvalue 阈值 pAdjustMethod = "BH", #Pvalue 矫正方法 
                readable = FALSE)        #TRUE 则展示SYMBOL，FALSE 则展示原来的ID（选false是因为不是所有gene都有symbol的)

ego_results <- as.data.frame(ego)                          ###生成的ego文件转换成data.frame格式即可。

write.table(ego_results, file = "F0_F4_ego.txt", quote = F)                    ###让保存的字符串不用“”引起来
pdf(file = "F0_F4_ego_BP.pdf")                                                                   ##打开一个PDF文件
barplot(ego, showCategory=20, x = "GeneRatio")                                ##把图画到这个PDF文件里
dotplot(ego)
x2 <- pairwise_termsim(ego)
emapplot(x2)               
cnetplot(ego, showCategory = 5)
dev.off()                                                                                                 ##关闭PDF

