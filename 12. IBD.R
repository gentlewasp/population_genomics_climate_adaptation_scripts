rm(list = ls()) #���������

####IBD_ALL###
library("PBSmapping")                ###

sites_ALL <- read.table("sample.coord.txt",header=T)[,2:3]  ###���뾭γ�����ݣ������ڵ�һ�У�γ���ڵڶ��С�
#attach(sites)                         ###��X,Y�̶�ΪR�����е�Ĭ�����ݡ�
colnames(sites_ALL) <- c("X","Y")         ###������ת������ϵ�����������꣨��������ϵ������γ�ȣ���ƽ�����꣨ͶӰ����ϵ����ͨ��ī���У�ͶӰ��ƽ��ֱ������ϵ��
attr(sites_ALL,"projection") <- "LL" 
attr(sites_ALL,"zone") <- NA
sites.utm_ALL <- convUL(sites_ALL)
Dgeo.dis_ALL <- dist(sites.utm_ALL)               ###���������֮��ľ���
#Dgeo.dis_ALL <- log(Dgeo.dis_ALL)
clima_ALL <- read.table("4bios.txt",header=T)  ###���뾭γ�����ݣ������ڵ�һ�У�γ���ڵڶ��С�
clim.pca_ALL <- prcomp(clima_ALL,scale. = T)

env.dis_ALL <- dist(clim.pca_ALL$x[,1:2])


library(tseries)
library(ade4)                      
FST_ALL <- read.matrix("FST(1-FST)_R.txt")  ###�����Ŵ����룬ȫ����
FST.dis_ALL <- as.dist(FST_ALL)                                     ###ȫ����ת���ɾ�����������ǣ�

                                                    ###IBD���

ibd_ALL <- mantel.randtest(FST.dis_ALL,Dgeo.dis_ALL,nrepet = 1000000)                        ###IBD������mantel test
ibd_ALL 


####IBD_field###
library("PBSmapping")                ###

sites_field <- read.table("sample.coord_field.txt",header=T)[,2:3]  ###���뾭γ�����ݣ������ڵ�һ�У�γ���ڵڶ��С�
#attach(sites)                         ###��X,Y�̶�ΪR�����е�Ĭ�����ݡ�
colnames(sites_field) <- c("X","Y")         ###������ת������ϵ�����������꣨��������ϵ������γ�ȣ���ƽ�����꣨ͶӰ����ϵ����ͨ��ī���У�ͶӰ��ƽ��ֱ������ϵ��
attr(sites_field,"projection") <- "LL" 
attr(sites_field,"zone") <- NA
sites.utm_field <- convUL(sites_field)
Dgeo.dis_field <- dist(sites.utm_field)               ###���������֮��ľ���
#Dgeo.dis_field <- log(Dgeo.dis_field)
clima_field <- read.table("4bios_field.txt",header=T)  ###���뾭γ�����ݣ������ڵ�һ�У�γ���ڵڶ��С�
clim.pca_field <- prcomp(clima_field,scale. = T)

env.dis_field <- dist(clim.pca_field$x[,1:2])

library(tseries)
library(ade4)                      
FST_field <- read.matrix("FST(1-FST)_R_field.txt")  ###�����Ŵ����룬ȫ����
FST.dis_field <- as.dist(FST_field)                                     ###ȫ����ת���ɾ�����������ǣ�

                                                 ###IBD���

ibd_field <- mantel.randtest(FST.dis_field,Dgeo.dis_field,nrepet = 1000000)                        ###IBD������mantel test
ibd_field

####IBD_greenhouse###
library("PBSmapping")                ###

sites_greenhouse <- read.table("sample.coord_greenhouse.txt",header=T)[,2:3]  ###���뾭γ�����ݣ������ڵ�һ�У�γ���ڵڶ��С�
#attach(sites)                         ###��X,Y�̶�ΪR�����е�Ĭ�����ݡ�
colnames(sites_greenhouse) <- c("X","Y")         ###������ת������ϵ�����������꣨��������ϵ������γ�ȣ���ƽ�����꣨ͶӰ����ϵ����ͨ��ī���У�ͶӰ��ƽ��ֱ������ϵ��
attr(sites_greenhouse,"projection") <- "LL" 
attr(sites_greenhouse,"zone") <- NA
sites.utm_greenhouse <- convUL(sites_greenhouse)
Dgeo.dis_greenhouse <- dist(sites.utm_greenhouse)               ###���������֮��ľ���
#Dgeo.dis_greenhouse <- log(Dgeo.dis_greenhouse)
clima_greenhouse <- read.table("4bios_greenhouse.txt",header=T)  ###���뾭γ�����ݣ������ڵ�һ�У�γ���ڵڶ��С�
clim.pca_greenhouse <- prcomp(clima_greenhouse,scale. = T)

env.dis_greenhouse <- dist(clim.pca_greenhouse$x[,1:2])

library(tseries)
library(ade4)                      
FST_greenhouse <- read.matrix("FST(1-FST)_R_greenhouse.txt")  ###�����Ŵ����룬ȫ����
FST.dis_greenhouse <- as.dist(FST_greenhouse)                                     ###ȫ����ת���ɾ�����������ǣ�



ibd_greenhouse <- mantel.randtest(FST.dis_greenhouse,Dgeo.dis_greenhouse,nrepet = 1000000)                        ###IBD������mantel test
ibd_greenhouse



library(ggplot2) #����ggplot2��
library(ggpmisc) #����ggpmisc��
library(dplyr)
l_gen_ALL <- as.vector(FST.dis_ALL)
l_geo_ALL <- as.vector(Dgeo.dis_ALL)
l_eco_ALL <- as.vector(env.dis_ALL)

l_gen_field <- as.vector(FST.dis_field)
l_geo_field <- as.vector(Dgeo.dis_field)
l_eco_field <- as.vector(env.dis_field)

l_gen_greenhouse <- as.vector(FST.dis_greenhouse)
l_geo_greenhouse <- as.vector(Dgeo.dis_greenhouse)
l_eco_greenhouse <- as.vector(env.dis_greenhouse)



#IBD��ͼ
pdf("0.5tp.IBD_R_fst_4bios_4.14.pdf",width = 4,height = 3)
ggplot() +

  geom_point(data=ibdibe_ALL,aes(x=l_geo_ALL, y=l_gen_ALL),color="black",alpha = 0.3,size=0.5) +
  geom_smooth(data=ibdibe_ALL,aes(x=l_geo_ALL, y=l_gen_ALL),method = "lm", se = TRUE,color="black",alpha = 0.2)+
  
  geom_point(data=ibdibe_field,aes(x=l_geo_field, y=l_gen_field),color="red",alpha = 0.3,size=0.5) +
  geom_smooth(data=ibdibe_field,aes(x=l_geo_field, y=l_gen_field),method = "lm", se = TRUE,color="red",alpha = 0.2)+
  
  geom_point(data=ibdibe_greenhouse,aes(x=l_geo_greenhouse, y=l_gen_greenhouse),color="blue",alpha = 0.3,size=0.5) +
  geom_smooth(data=ibdibe_greenhouse,aes(x=l_geo_greenhouse, y=l_gen_greenhouse),method = "lm", se = TRUE,color="blue",alpha = 0.2)+
  

  theme_bw()+theme(panel.grid = element_blank())
dev.off()


























