---
title: Spinal cord injury changes the structure and functional potential of gut bacterial
  and viral communities in a spinal-level dependent manner
author: "Jingjie Du"
date: "4/8/2020"
output:
  html_document:
    self_contained: yes 
    highlight: kate
    theme: yeti
    toc: yes
  pdf_document:
    keep_tex: yes
    toc: yes
  word_document:
    toc: yes
---

```{r setup, include=FALSE}
## knitr options - don't delete!
knitr::opts_chunk$set(
  fig.path = "Figures/"
)
```


```{r}
library(pracma)
library(vegan)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library("affy") 
library("limma")
library("multtest")
library("DESeq2")
library(scales)
library(tidyr)
library("hgu95av2cdf")
library("pheatmap")
library("dendextend")
library(funrar)
library(edgeR)

setwd("~/Desktop/spinal cord injury project/Ion torrent paper")

Colvec_expedition=read.csv("Color_vector.csv",header=T,row.names=1,sep=",")
```

#Figure 1. Intestinal bacterial community composition was disturbed after spinal cord injury

## B. Principal coordinate analysis (PCoA) of Bray-Curtis distances using single copy gene L2
```{r}
L2=read.csv("otu_table.L2.csv",header=T,row.names=1,sep=",")
L2=t(L2) 
L2=t(L2) 
#calculate the relative abundances of bacterial strains
L2_abundances=make_relative(L2)

colvec_all <- c("#999999", "#999999","#999999","#999999","#999999","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00", "#56B4E9","#56B4E9","#56B4E9","#56B4E9","#56B4E9")
# Applying a log transformation
d.clr_L2= log2(L2_abundances+1)
# Calculating Bray-Curtis dissimilarity
Bray_L2=vegdist(d.clr_L2,"bray")
# ANOSIM statistics
L2_anosim=anosim(Bray_L2, permutations = 999, grouping=Colvec_expedition$Zones)
summary(L2_anosim)
#ANOSIM statistic R: 0.3973 , Significance: 0.001 

# Applying a non-constrained ordination (PCoA) on the dissimilarity matrix 
L2_bray=capscale(Bray_L2~-1)
# Extracting the percentage explained by the first two dimensions and automatically adding them to the axes titles
L2_bray_eig = eigenvals(L2_bray)
percentage_variance_explained <- L2_bray_eig / sum(L2_bray_eig)
sum_percentage_variance_explained <- cumsum(L2_bray_eig / sum(L2_bray_eig))
xlabel= as.numeric(format(round((percentage_variance_explained[1]*100), 2), nsmall = 2))
xlabel= sprintf("%.2f %%", xlabel)
xlabel= paste ("PCo1 (", xlabel, "of variation explained", ")")
ylabel= as.numeric(format(round((percentage_variance_explained[2]*100), 2), nsmall = 2))
ylabel= sprintf("%.2f %%", ylabel)
ylabel= paste ("PCo2 (", ylabel,"of variation explained", ")")

# Plotting the figure
# Adding the axes, grid, and other aestethics
plot(L2_bray, type="n", xlab="", xlim=c(-1,1),ylab="", tck = -0.01, mgp = c(3, 0.2, 0), 
     xaxp  = c(-4, 4, 8), panel.first=grid(col = "white",lty=0))
title(ylab=ylabel, line=2, cex.lab=1)
title(xlab=xlabel, line=2, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
points(L2_bray,cex= 2.5,pch=21, col="black", bg= colvec_all, lwd = 1)
# Adding a legend
legend(-2,1, pt.cex=2.5 , pt.lwd = 1, c("Lam", "T10", "T4"), bty = "n", pch = 21, col="black", pt.bg = c("#999999","#E69F00","#56B4E9"), cex = 1.5)
ordihull(L2_bray,Colvec_expedition$Zones, display = "sites", col = NULL,label=T)

```
##C.Bray_Curtis dissimalarity between the control group and SCI groups
```{r}
Bray_L2_SCI=read.csv("Bray_L2_SCI.csv",header=T,sep=",")

pairwise.wilcox.test(Bray_L2_SCI$Bray_curtis, Bray_L2_SCI$Zones,p.adjust.method = "fdr")
#p-value = 0.0052
boxplot_Bray_L2_SCI=ggplot(Bray_L2_SCI,  aes(x=Zones, y=Bray_curtis)) + geom_boxplot(fill=c( "#E69F00", "#56B4E9"))+geom_jitter(width = 0.06)+labs(title=NULL,x=NULL,y =c( "Bray_curtis dissimilarity"))+theme_classic()+geom_point(shape=16,fill="black",size=2)

boxplot_Bray_L2_SCI1=boxplot_Bray_L2_SCI+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20,face = "bold"))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
#Including boxplots
p=boxplot_Bray_L2_SCI1+scale_y_continuous(expand = c(0,0), limits=c(0.3,0.9))+geom_segment(y=0.83,yend=0.83,x=1,xend=2)+geom_text(x=1.5,y=0.83,label="**",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
p

```


```{r}
#Read bacterial abundances in phylum level
phylum_singleM=read.csv("15_ion_annatable_phylum_singleM.csv",header=T,row.names=1,sep=",") 
phylum_singleM1=t(phylum_singleM)
phylum_singleM2=as.data.frame(phylum_singleM1)
phylum_singleM2$Sample=factor(rownames(phylum_singleM2))
phylum_singleM2_m=melt(phylum_singleM2)

rare_phylum_singleM=data.frame(phylum_singleM2$` p__Actinobacteria`,phylum_singleM2$` p__Proteobacteria`,phylum_singleM2$` p__Firmicutes_B`,phylum_singleM2$others)
rownames(rare_phylum_singleM)=rownames(phylum_singleM2)
rare_phylum_singleM$Sample=factor(rownames(rare_phylum_singleM))
rare_phylum_singleM_m=melt(rare_phylum_singleM)

```

## D.Intestinal bacterial community composition at the phylum level.
```{r}
phylum_singleM_p=ggplot(phylum_singleM2_m)+geom_bar(aes(x=Sample,y=value,fill=variable),stat="identity")+labs(title = NULL,x=NULL,y = "Relative abundance")+ scale_y_continuous(limits = c(0,1.00), expand = c(0, 0))+theme_bw()

# legend labels
phylum_singleM_p=phylum_singleM_p+geom_segment(y=0,yend=1.00,x=5.5,xend=5.5)+geom_segment(y=0,yend=1.00,x=10.5,xend=10.5)

phylum_singleM_p+theme(plot.title = element_text(size = 20,hjust=0.5))+scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))+theme(axis.text.x = element_text(color="black",size = 12))+guides(fill=guide_legend(title=NULL))+theme(legend.title = element_text(size=18))+ theme(legend.text = element_text(size=14),legend.position="right")+guides(fill=guide_legend("Phylum"))+scale_fill_brewer(palette = "Set2")
#On the x-axis, the numbers 1–15 represent individual mice within each group. numbers 1-5 represents Lam, numbers 6-10 represents T10 and numbers 11-15 represents T4
```

```{r}
rare_phylum_singleM_p=ggplot(rare_phylum_singleM_m)+geom_bar(aes(x=Sample,y=value,fill=variable),stat="identity")+labs(title = NULL,x=NULL,y = "Relative abundance")+theme_bw()
# legend labels
#phylum_singleM_p=phylum_singleM_p+geom_segment(y=0,yend=1.00,x=5.5,xend=5.5)+geom_segment(y=0,yend=0.1,x=10.5,xend=10.5)

rare_phylum_singleM_p+theme(plot.title = element_text(size = 20,hjust=0.5))+scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))+theme(axis.text.x = element_text(color="black",size = 12))+guides(fill=guide_legend(title=NULL))+theme(legend.title = element_text(size=18))+ theme(legend.text = element_text(size=14),legend.position="right")+guides(fill=guide_legend("Phylum"))+scale_colour_manual(values = c("#A6D854","#FFD92F","#E5C494","#B3B3B3"))
#On the x-axis, the numbers 1–15 represent individual mice within each group. numbers 1-5 represents Lam, numbers 6-10 represents T10 and numbers 11-15 represents T4
```

```{r}
##Wilcoxon rank sum test for bacterial phyla
clin.data <- data.frame(ID = 1:ncol(phylum_singleM), group = rep(c("Lam", "T10","T4"), each=5))
clin.data

idx_Lam<-sample(which(clin.data $group=="Lam"), 5, replace=FALSE) 
idx_T10 <- sample(which(clin.data $group=="T10"), 5, replace=FALSE) 
idx_T4 <- sample(which(clin.data $group=="T4"), 5, replace=FALSE)

#Wilcoxon rank sum test between Lam and T4 at phylum level
phylum_data_test <- phylum_singleM[,c(idx_Lam, idx_T4)]
clin.data_phylum_sub <-  data.frame(ID = 1:ncol(phylum_data_test), group = rep(c("Lam","T4"), each=5))
phylum_data_test1=t(phylum_data_test)
phylum_data_test2=cbind(phylum_data_test1,clin.data_phylum_sub)

wilcox.test(phylum_data_test2$` p__Bacteroidetes`~phylum_data_test2$group)
#p-value = 0.22
wilcox.test(phylum_data_test2$` p__Firmicutes`~phylum_data_test2$group)
#p-value = 0.8413
wilcox.test(phylum_data_test2$` p__Firmicutes_A`~phylum_data_test2$group)
#p-value = 0.1508
wilcox.test(phylum_data_test2$` p__Actinobacteria`~phylum_data_test2$group)
#p-value = 0.04653 *
wilcox.test(phylum_data_test2$` p__Proteobacteria`~phylum_data_test2$group)
#p-value = 1

##Wilcoxon rank sum test obetween Lam and T10 at phylum level
phylum_data_test4 <- phylum_singleM[,c(idx_Lam, idx_T10)]
clin.data_phylum_sub1 <-  data.frame(ID = 1:ncol(phylum_data_test), group = rep(c("Lam","T10"), each=5))
phylum_data_test5=t(phylum_data_test4)
phylum_data_test6=cbind(phylum_data_test5,clin.data_phylum_sub1)

wilcox.test(phylum_data_test6$` p__Bacteroidetes`~phylum_data_test6$group)
#p-value = 1
wilcox.test(phylum_data_test6$` p__Firmicutes`~phylum_data_test6$group)
#p-value = 0.03175 *
wilcox.test(phylum_data_test6$` p__Firmicutes_A`~phylum_data_test6$group)
#p-value = 0.8413
wilcox.test(phylum_data_test6$` p__Actinobacteria`~phylum_data_test6$group)
#p-value = 0.0278 *
wilcox.test(phylum_data_test6$` p__Proteobacteria`~phylum_data_test6$group)
# p-value = 0.08326
```

```{r}
#Wilcoxon rank sum test corrected by BH
singleM_phylum=cbind(phylum_singleM1,Colvec_expedition)
colnames(singleM_phylum)=c("Bacteroidetes","Firmicutes","Firmicutes_A","Actinobacteria","Proteobacteria","Firmicutes_B","others","Zones")
pairwise.wilcox.test(singleM_phylum$Actinobacteria, singleM_phylum$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_phylum$Firmicutes, singleM_phylum$Zones,p.adjust.method = "fdr")
```

##E.Boxplots showing the relative abundance of the phyla Firmicutes 
```{r}
#Boxplots for Firmicutes
Firmicutes=ggplot(singleM_phylum, aes(x=Zones, y=Firmicutes)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+labs(title = NULL,x=NULL,y = "Firmicutes \n relative abundances")+theme_classic()
Firmicutes_boxplot=Firmicutes+theme(plot.title = element_text(size = 20,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 16),axis.title.y=element_text(size=18))
Firmicutes_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,0.25))+geom_segment(y=0.23,yend=0.23,x=1,xend=2)+geom_text(x=1.5,y=0.23,label="*",size=12)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

##F.Boxplots showing the relative abundance of the phyla Actinobacteria
```{r}
Actinobacteria=ggplot(singleM_phylum, aes(x=Zones, y=Actinobacteria)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+labs(title = NULL,x=NULL,y = "Actinobacteria \n relative abundances")+theme_classic()
Actinobacteria_boxplot=Actinobacteria+theme(plot.title = element_text(size = 20,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 16),axis.title.y=element_text(size=18))
Actinobacteria_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,0.03))+geom_segment(y=0.025,yend=0.025,x=1,xend=2)+geom_text(x=1.5,y=0.025,label="*",size=12)+geom_segment(y=0.027,yend=0.027,x=1,xend=3)+geom_text(x=2,y=0.027,label="*",size=12)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
#*p = 0.05 by Wilcoxon rank sum test . Each dot in the boxplot represents an individual mouse sample, and each group (Lam, T4, T10) contains five samples. 
```

```{r}
#Read bacterial abundances in genus
genus_singleM=read.csv("15_ion_annatable_genus_singleM.csv",header=T,row.names=1,sep=",") 

##Wilcoxon rank sum test 

#Lam and T4
genus_singleM1<- genus_singleM[,c(idx_Lam, idx_T4)]
clin.sub_genus_singleM1 <-  data.frame(ID = 1:ncol(genus_singleM1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_genus_singleM1$group)
  return(res$p.value)
}
genus_singleM.wilcox1 <-apply(as.matrix(genus_singleM1), 1, fx)
genus_singleM.wilcox1 <- data.frame(genus_singleM.wilcox1 ) 
colnames(genus_singleM.wilcox1) <- c("p.value")
write.csv(genus_singleM.wilcox1,file="genus_wilcox0.05_Lam_T4.csv")
rownames(genus_singleM.wilcox1)[genus_singleM.wilcox1$p.value< 0.01]
#**,p<0.01, "Turicibacter" ; "Lactobacillus"
rownames(genus_singleM.wilcox1)[genus_singleM.wilcox1$p.value< 0.05]
#*,p<0.05, "Ruminococcus_A"; "Bacteroides" ;"Neglecta"; "UBA9475";"Eubacterium_F";" CAG-1031"     

#Lam and T10
genus_singleM2 <- genus_singleM[,c(idx_Lam, idx_T10)]
clin.sub_genus_singleM2 <-  data.frame(ID = 1:ncol(genus_singleM2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_genus_singleM2$group)
  return(res$p.value)
}
genus_singleM.wilcox2 <-apply(as.matrix(genus_singleM2), 1, fx)
genus_singleM.wilcox2 <- data.frame(genus_singleM.wilcox2) 
colnames(genus_singleM.wilcox2) <- c("p.value")
write.csv(genus_singleM.wilcox2,file="genus_wilcox0.05_Lam_T10.csv")
rownames(genus_singleM.wilcox2)[genus_singleM.wilcox2$p.value< 0.01]
#**, p<0.01, "Eubacterium_R"
rownames(genus_singleM.wilcox2)[genus_singleM.wilcox2$p.value< 0.05]
#*, p<0.05, "Slackia";"Pseudooceanicola";"Ruminococcus";"Flavonifractor";"Lachnospira";"Clostridium_M";"Turicibacter";"Lactobacillus"

####T4 and T10
genus_singleM3 <- genus_singleM[,c(idx_T4, idx_T10)]
clin.sub_genus_singleM3 <-  data.frame(ID = 1:ncol(genus_singleM3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_genus_singleM3$group)
  return(res$p.value)
}
genus_singleM.wilcox3 <-apply(as.matrix(genus_singleM3), 1, fx)
genus_singleM.wilcox3 <- data.frame(genus_singleM.wilcox3) 
colnames(genus_singleM.wilcox3) <- c("p.value")
write.csv(genus_singleM.wilcox3,file="genus_wilcox0.05_T4_T10.csv")
rownames(genus_singleM.wilcox3)[genus_singleM.wilcox3$p.value< 0.01]
#**,p<0.01,  "CAG-791","Eubacterium_I"
rownames(genus_singleM.wilcox3)[genus_singleM.wilcox3$p.value< 0.05]
#*,p<0.05, "Weissella","Clostridium_M"
```

```{r}
#pairwise.wilcox.test with BH correction
genus_singleM4=genus_singleM*100
genus_singleM5=t(genus_singleM4)
singleM_genus=cbind(genus_singleM5,Colvec_expedition)

pairwise.wilcox.test(singleM_genus$` CAG-1031`, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Lactobacillus, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Turicibacter, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Bacteroides, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Weissella, singleM_genus$Zones,p.adjust.method = "fdr")
#members in class Clostridia
pairwise.wilcox.test(singleM_genus$Eubacterium_R, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Eubacterium_F, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$UBA9475, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Neglecta, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Lachnospira, singleM_genus$Zones,p.adjust.method = "fdr")
#rare genera
pairwise.wilcox.test(singleM_genus$Clostridium_M, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Eubacterium_I, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$`CAG-791`, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Slackia, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Pseudooceanicola, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Ruminococcus, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Flavonifractor, singleM_genus$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(singleM_genus$Ruminococcus_A, singleM_genus$Zones,p.adjust.method = "fdr")
```



#Figure 2. Genus-level bacterial abundances are altered after SCI.

##A.Differential abundances of bacterial genera (p<0.05 by Wilcoxon signed-rank test,FDR<0.05)
```{r  fig.height = 12, fig.width = 5, fig.align = "center"}
#B.Differential abundances of bacterial genera (p<0.05 by Wilcoxon signed-rank test) in either two groups are indicated in red. Each row representing a unique bacterial genus was Z-score normalized. Bacterial genera on the y-axis are clustered using Euclidean distances.  
#use pheatmap

genus_singleM_FDR_0.05=data_frame(singleM_genus$Clostridium_M,singleM_genus$Eubacterium_I,singleM_genus$Slackia,singleM_genus$`CAG-791`,singleM_genus$Flavonifractor,singleM_genus$Lactobacillus,singleM_genus$Turicibacter,singleM_genus$` CAG-1031`,singleM_genus$Bacteroides,singleM_genus$Eubacterium_R,singleM_genus$Eubacterium_F)
colnames(genus_singleM_FDR_0.05)=c("Clostridium_M","Eubacterium_I","Slackia","CAG-791","Flavonifractor","Lactobacillus","Turicibacter","CAG-1031","Bacteroides","Eubacterium_R","Eubacterium_F")
rownames(genus_singleM_FDR_0.05)=rownames(singleM_genus)

genus_singleM_FDR_0.05_t=t(genus_singleM_FDR_0.05)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

genus_singleM_FDR_0.05_norm <- t(apply(genus_singleM_FDR_0.05_t, 1, cal_z_score))

pheatmap(genus_singleM_FDR_0.05_norm,cluster_cols=FALSE, legend = TRUE,border_color = NA)

```


##B-D. Relative abundances of select groups indicating that  CAG-1031, Lactobacillus, and Turicibacter decreased after SCI 
```{r} 
#Box-plots for significant genus

#A. CAG-1031
CAG=ggplot(singleM_genus, aes(x=Zones, y=` CAG-1031`)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "CAG-1031 \n relative abundances(%)")+theme_classic()
CAG_boxplot=CAG+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))

CAG_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,40))+geom_segment(y=36,yend=36,x=1,xend=3)+geom_text(x=2,y=36,label="*",size=12)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

#B. Lactobacillus
Lactobacillus=ggplot(singleM_genus, aes(x=Zones, y=Lactobacillus)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Lactobacillus \n relative abundances(%)")+theme_classic()

Lactobacillus_boxplot=Lactobacillus+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Lactobacillus_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,10))+geom_segment(y=8.8,yend=8.8,x=1,xend=2)+geom_segment(y=9.5,yend=9.5,x=1,xend=3)+geom_text(x=1.5,y=8.8,label="*",size=12)+geom_text(x=2,y=9.5,label="**",size=12)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

#Turicibacter
Turicibacter=ggplot(singleM_genus, aes(x=Zones, y=Turicibacter)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Turicibacter \n relative abundances(%)")+theme_classic()

Turicibacter_boxplot=Turicibacter+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Turicibacter_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,4.6))+geom_segment(y=3.8,yend=3.8,x=1,xend=2)+geom_segment(y=4.2,yend=4.2,x=1,xend=3)+geom_text(x=1.5,y=3.8,label="*",size=12)+geom_text(x=2,y=4.2,label="**",size=12)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

##E-G. Bacteroides, Eubacterium_F and Eubacterium_R increased after SCI  compared to Lam controls. 
```{r}
#e. Bacteroides
Bacteroides=ggplot(singleM_genus, aes(x=Zones, y=Bacteroides)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Bacteroides \n relative abundances(%)")+theme_classic()

Bacteroides_boxplot=Bacteroides+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Bacteroides_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,7))+geom_segment(y=6.4,yend=6.4,x=1,xend=3)+geom_text(x=2,y=6.4,label="*",size=12)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

#f. Eubacterium_F
Eubacterium_F=ggplot(singleM_genus, aes(x=Zones, y=Eubacterium_F)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Eubacterium_F \n relative abundances(%)")+theme_classic()

Eubacterium_F_boxplot=Eubacterium_F+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Eubacterium_F_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,5))+geom_segment(y=4.5,yend=4.5,x=1,xend=3)+geom_text(x=2,y=4.5,label="*",size=12)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))


#g.Eubacterium_R
Eubacterium_R=ggplot(singleM_genus, aes(x=Zones, y=Eubacterium_R)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Eubacterium_R \n relative abundances(%)")+theme_classic()
Eubacterium_R_boxplot=Eubacterium_R+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Eubacterium_R_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,4))+geom_segment(y=3.6,yend=3.6,x=1,xend=2)+geom_text(x=1.5,y=3.6,label="**",size=12)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

```


#Figure 3. Species-level bacterial abundances are altered after SCI. 
```{r}
#read MAGs 
mOTUs=read.table("microbial_bins_60percent_complete.txt",header=T,row.names=1,sep="\t")
mOTUs$Maxbins=as.factor(rownames((mOTUs)))
mOTUs_species=read.table("MAGs_species.txt",header=T,row.names=1,sep="\t")
mOTUs_species$Maxbins=as.factor(rownames((mOTUs_species)))
mOTUs_species1=merge(mOTUs_species,mOTUs,by="Maxbins")
mOTUs_species2=mOTUs_species1[,-1]
row.names(mOTUs_species2)=mOTUs_species2[,1]
mOTUs_species3=mOTUs_species2[,-1]

#Lam and T4
MAGs_species1 <- mOTUs_species3[,c(idx_Lam, idx_T4)]
clin.sub_MAGs_species1 <-  data.frame(ID = 1:ncol(MAGs_species1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_MAGs_species1$group)
  return(res$p.value)
}

MAGs_species.wilcox1 <-apply(as.matrix(MAGs_species1), 1, fx)
MAGs_species.wilcox1 <- data.frame(MAGs_species.wilcox1 ) 
colnames(MAGs_species.wilcox1) <- c("p.value")
write.csv(MAGs_species.wilcox1, file = "MAGs_species.wilcox0.05_Lam_T4.csv")
rownames(MAGs_species.wilcox1)[MAGs_species.wilcox1$p.value< 0.01]
#**, "Lactococcus lactis_A"    "Lactobacillus johnsonii"
rownames(MAGs_species.wilcox1)[MAGs_species.wilcox1$p.value< 0.05]
#*,"Weissella cibaria"

#Lam and T10
MAGs_species2 <- mOTUs_species3[,c(idx_Lam, idx_T10)]
clin.sub_MAGs_species2 <-  data.frame(ID = 1:ncol(MAGs_species2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_MAGs_species2$group)
  return(res$p.value)
}

MAGs_species.wilcox2 <-apply(as.matrix(MAGs_species2), 1, fx)
MAGs_species.wilcox2 <- data.frame(MAGs_species.wilcox2 ) 
colnames(MAGs_species.wilcox2) <- c("p.value")
write.csv(MAGs_species.wilcox2, file = "MAGs_species.wilcox0.05_Lam_T10.csv")
rownames(MAGs_species.wilcox2)[MAGs_species.wilcox2$p.value< 0.01]
#**,"Weissella cibaria"
rownames(MAGs_species.wilcox2)[MAGs_species.wilcox2$p.value< 0.05]
#*,"Lactobacillus johnsonii"

#T10 and T4
MAGs_species3 <- mOTUs_species3[,c(idx_T4, idx_T10)]
clin.sub_MAGs_species3 <-  data.frame(ID = 1:ncol(MAGs_species3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_MAGs_species3$group)
  return(res$p.value)
}

MAGs_species.wilcox3 <-apply(as.matrix(MAGs_species3), 1, fx)
MAGs_species.wilcox3 <- data.frame(MAGs_species.wilcox3 ) 
colnames(MAGs_species.wilcox3) <- c("p.value")
write.csv(MAGs_species.wilcox3, file = "MAGs_species.wilcox0.05_T4_T10.csv")
rownames(MAGs_species.wilcox3)[MAGs_species.wilcox3$p.value< 0.05]
#**,"Lactococcus lactis_A"
rownames(MAGs_species.wilcox3)[MAGs_species.wilcox3$p.value< 0.05]
#*, "Weissella cibaria" 
```

```{r}
mOTUs_species4=t(mOTUs_species3)
mOTUs_species5=cbind(mOTUs_species4, Colvec_expedition)

#Wilcoxon Rank Sum test with 'fdr'('BH') correction for an FDR adjusted p-value (or q-value)
pairwise.wilcox.test(mOTUs_species5$`Lactobacillus johnsonii`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.009`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.121`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`Weissella cibaria`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`Lactococcus lactis_A`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`Bacteroides thetaiotaomicron`, mOTUs_species5$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(mOTUs_species5$`maxbin.004`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.007`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.011`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.014`, mOTUs_species5$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(mOTUs_species5$`maxbin.017`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.019`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.020`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.026`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.029`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.033`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.034`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.043`, mOTUs_species5$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(mOTUs_species5$`maxbin.044`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.056`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.057`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.062`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.067`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.071`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.087`, mOTUs_species5$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(mOTUs_species5$`maxbin.123`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.124`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.133`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.135`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.141`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.151`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.152`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.157`, mOTUs_species5$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(mOTUs_species5$`maxbin.166`, mOTUs_species5$Zones,p.adjust.method = "fdr")

```

```{r}
#differential MAGs
MAGS_differential=data_frame(mOTUs_species5$`Lactobacillus johnsonii`,mOTUs_species5$`maxbin.009`,mOTUs_species5$`maxbin.121`, mOTUs_species5$`Weissella cibaria`, mOTUs_species5$`Lactococcus lactis_A`,mOTUs_species5$`maxbin.004`,mOTUs_species5$`maxbin.007`,mOTUs_species5$`maxbin.011`,mOTUs_species5$`maxbin.014`,mOTUs_species5$`maxbin.017`,mOTUs_species5$`maxbin.019`,mOTUs_species5$`maxbin.020`,mOTUs_species5$`maxbin.026`,mOTUs_species5$`maxbin.029`,mOTUs_species5$`maxbin.033`,mOTUs_species5$`maxbin.034`,mOTUs_species5$`maxbin.043`,mOTUs_species5$`maxbin.044`,mOTUs_species5$`maxbin.056`,mOTUs_species5$`maxbin.057`,mOTUs_species5$`maxbin.062`,mOTUs_species5$`maxbin.067`,mOTUs_species5$`maxbin.071`,mOTUs_species5$`maxbin.087`,mOTUs_species5$`maxbin.123`,mOTUs_species5$`maxbin.124`,mOTUs_species5$`maxbin.133`,mOTUs_species5$`maxbin.135`,mOTUs_species5$`maxbin.141`,mOTUs_species5$`maxbin.151`,mOTUs_species5$`maxbin.152`,mOTUs_species5$`maxbin.157`,mOTUs_species5$`maxbin.166`)

rownames(MAGS_differential)=rownames(mOTUs_species5)

MAGS_differential_t=t(MAGS_differential)
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(MAGS_differential_t)
pheatmap(MAGS_differential_t,cluster_rows=T,cluster_cols=FALSE,annotation_col = my_sample_col)

```


##A.Differential MAGs with FDR<0.05 
```{r , fig.height = 12, fig.width = 10, fig.align = "center"}
#Only keep MAGs with FDR<0.05 

mOTUs_species_new=read.table("MAGs_species_FDR_0.05.txt",header=F,row.names=1,sep="\t")
mOTUs_species_new$Maxbins=as.factor(rownames((mOTUs_species_new)))
mOTUs_species_new1=merge(mOTUs_species_new,mOTUs,by="Maxbins")
mOTUs_species_new2=mOTUs_species_new1[,-1]
row.names(mOTUs_species_new2)=mOTUs_species_new2[,1]
mOTUs_species_new3=mOTUs_species_new2[,-1]

mOTUs_species_new_norm <- t(apply(mOTUs_species_new3, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(mOTUs_species_new3)
#pdf(file="Figure S3_new",height = 12,width = 10)
pheatmap(mOTUs_species_new_norm,cluster_cols=FALSE, annotation_col = my_sample_col,legend = TRUE,cluster_rows = T,border_color = NA,cellwidth = 12, cellheight = 15)

```


##B-D, Beneficial bacterial species(MAGs) decreased after SCI
###B. Lactobacillus johnsonii
```{r}
Lactobacillus_johnsonii=ggplot(mOTUs_species5, aes(x=Zones, y=`Lactobacillus johnsonii`)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "Lactobacillus johnsonii",x=NULL,y ="RPKM")+theme_classic()
Lactobacillus_johnsonii_boxplot=Lactobacillus_johnsonii+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 20),axis.title.y=element_text(size=20))
Lactobacillus_johnsonii_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,45))+geom_segment(y=40,yend=40,x=1,xend=2)+geom_segment(y=42,yend=42,x=1,xend=3)+geom_text(x=2,y=42,label="**",size=10)+geom_text(x=1.5,y=40,label="*",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

###C. maxbin.009 CAG-1031
```{r}
maxbin009=ggplot(mOTUs_species5, aes(x=Zones, y=`maxbin.009`)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+labs(title = "CAG-1031 maxbin.009",x=NULL,y ="RPKM")+theme_classic()

maxbin009_boxplot=maxbin009+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 20),axis.title.y=element_text(size=20))
maxbin009_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,82))+geom_segment(y=78,yend=78,x=1,xend=3)+geom_text(x=2,y=78,label="**",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

###D. maxbin.121 CAG-1031
```{r}
maxbin121=ggplot(mOTUs_species5, aes(x=Zones, y=`maxbin.121`)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+labs(title = "CAG-1031 maxbin.121",x=NULL,y ="RPKM")+theme_classic()
maxbin121_boxplot=maxbin121+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 20),axis.title.y=element_text(size=20))
maxbin121_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,64))+geom_segment(y=60,yend=60,x=1,xend=3)+geom_text(x=2,y=60,label="**",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

##E-G, Potential pathogenic bacterial species(MAGs) increased after SCI

###E. Weissella cibaria
```{r}
Weissella_cibaria=ggplot(mOTUs_species5, aes(x=Zones, y=`Weissella cibaria`)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+labs(title = "Weissella cibaria",x=NULL,y ="RPKM")+theme_classic()
Weissella_cibaria_boxplot=Weissella_cibaria+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 20),axis.title.y=element_text(size=20))
Weissella_cibaria_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,41))+geom_segment(y=37,yend=37,x=2,xend=3)+geom_segment(y=39,yend=39,x=1,xend=3)+geom_text(x=2,y=39,label="*",size=12)+geom_text(x=2.5,y=37,label="*",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

###F. Lactococcus lactis_A
```{r}
Lactococcus_lactis_A=ggplot(mOTUs_species5, aes(x=Zones, y=`Lactococcus lactis_A`)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+labs(title = "Lactococcus lactis_A",x=NULL,y ="RPKM")+theme_classic()
Lactococcus_lactis_A_boxplot=Lactococcus_lactis_A+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 20),axis.title.y=element_text(size=20))
  
  Lactococcus_lactis_A_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,40))+geom_segment(y=36,yend=36,x=2,xend=3)+geom_segment(y=38,yend=38,x=1,xend=3)+geom_text(x=2,y=38,label="**",size=10)+geom_text(x=2.5,y=36,label="**",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

###G. Bacteroides thetaiotaomicron
```{r}
Bacteroides_thetaiotaomicron=ggplot(mOTUs_species5, aes(x=Zones, y=`Bacteroides thetaiotaomicron`)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "Bacteroides thetaiotaomicron",x=NULL,y ="RPKM")+theme_classic()
Bacteroides_thetaiotaomicron_boxplot=Bacteroides_thetaiotaomicron+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 20),axis.title.y=element_text(size=20))
Bacteroides_thetaiotaomicron_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,22))+geom_segment(y=19,yend=19,x=1,xend=3)+geom_text(x=2,y=20,label="p=0.055",size=6)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

#Figure 4. Predicted metabolic pathways are different between healthy and spinal cord injury animals. 
##A.  predicted protein clusters are different between the control group Lam and SCI disease groups.
```{r}
unique_PCs=read.csv("PCs_RPKM.csv",header=T,sep=",")
rownames(unique_PCs)=unique_PCs[,1]
unique_PCs=unique_PCs[,-1]
t_unique_PCs=t(unique_PCs)
# Applying a cube root transformation
d.clr_PCs = log2(t_unique_PCs+1)
# Calculating Bray-Curtis dissimilarity
Bray_PCs=vegdist(d.clr_PCs,"bray")
# Applying a non-constrained ordination (PCoA) on the dissimilarity matrix 
unique_PCs_bray=capscale(Bray_PCs~-1)
# Extracting the percentage explained by the first two dimensions and automatically adding them to the axes titles
unique_PCs_bray_eig = eigenvals(unique_PCs_bray)
percentage_variance_explained <- unique_PCs_bray_eig/sum(unique_PCs_bray_eig)
sum_percentage_variance_explained <- cumsum(unique_PCs_bray_eig / sum(unique_PCs_bray_eig))
xlabel= as.numeric(format(round((percentage_variance_explained[1]*100), 2), nsmall = 2))
xlabel= sprintf("%.2f %%", xlabel)
xlabel= paste ("PCo1 (", xlabel, "of variation explained", ")")
ylabel= as.numeric(format(round((percentage_variance_explained[2]*100), 2), nsmall = 2))
ylabel= sprintf("%.2f %%", ylabel)
ylabel= paste ("PCo2 (", ylabel,"of variation explained", ")")

plot(unique_PCs_bray, type="n", xlab="", ylab="",ylim=c(-1, 1),xlim=c(-1,1),cex.axis=0.6, tck = -0.01, mgp = c(3, 1, 0), xaxp  = c(-2, 2, 4), panel.first=grid(col = "white",lty=0))

title(ylab=ylabel, line=2, cex.lab=1)
title(xlab=xlabel, line=2, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)

colvec_all <- c("#999999", "#999999","#999999","#999999","#999999","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00", "#56B4E9","#56B4E9","#56B4E9","#56B4E9","#56B4E9")

points(unique_PCs_bray,cex= 2.5,pch=21, col="black", bg= colvec_all, lwd = 1)
# Adding a legend
legend(-1.5,1, pt.cex=2.5 , pt.lwd = 1,c("Lam", "T10", "T4"),bty = "n" ,pch = 21,col="black",pt.bg = c("#999999", "#E69F00", "#56B4E9"), cex = 1.5)

ordihull(unique_PCs_bray, Colvec_expedition$Zones, display = "sites", label = T)

# ANOSIM statistics
PCs_anosim=anosim(Bray_PCs, permutations = 999, grouping=Colvec_expedition$Zones)
summary(PCs_anosim)
#ANOSIM statistic R: 0.4151 
#Significance: 0.001 
# Each data point indicates an individual mouse sample.
```

##B.Between_group Bray_curtis dissimilarity
```{r}

Bray_PCs_SCI=read.csv("Bray_PCs_SCI.csv",header=T,sep=",")
colnames(Bray_PCs_SCI)=c("Bray_curtis","Zones")
wilcox.test(Bray_PCs_SCI$Bray_curtis~Bray_PCs_SCI$Zone)
pairwise.wilcox.test(Bray_PCs_SCI$Bray_curtis, Bray_PCs_SCI$Zones,p.adjust.method = "fdr")
#p-value = 0.043

boxplot_Bray_PCs_SCI=ggplot(Bray_PCs_SCI,  aes(x=Zones, y=Bray_curtis)) + geom_boxplot(fill=c( "#E69F00", "#56B4E9"))+geom_jitter(width = 0.06)+labs(title=NULL,x=NULL,y =c( "Bray_curtis dissimilarity"))+theme_classic()+geom_point(shape=16,fill="black",size=2)

boxplot_Bray_PCs_SCI1=boxplot_Bray_PCs_SCI+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20,face = "bold"))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
#Including boxplots
boxplot_Bray_PCs_SCI1+scale_y_continuous(expand = c(0,0), limits=c(0.2,0.4))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +geom_segment(y=0.38,yend=0.38,x=1,xend=2)+geom_text(x=1.5,y=0.38,label="*",size=10)

```




```{r}
selected_functions=read.csv("selected_functions_MAGs.csv",header=T,row.names=NULL,sep=",")
row.names(selected_functions)=selected_functions[,1]
selected_functions=selected_functions[,-1]

#Lam and T4
selected_functions1 <- selected_functions[,c(idx_Lam, idx_T4)]
clin.sub_selected_functions1 <-  data.frame(ID = 1:ncol(selected_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_selected_functions1$group)
  return(res$p.value)
}

selected_functions.wilcox1 <-apply(as.matrix(selected_functions1), 1, fx)
selected_functions.wilcox1 <- data.frame(selected_functions.wilcox1 ) 
colnames(selected_functions.wilcox1) <- c("p.value")
write.csv(selected_functions.wilcox1, file = "selected_functions.wilcox0.05_Lam_T4.csv")
rownames(selected_functions.wilcox1)[selected_functions.wilcox1$p.value< 0.01]
rownames(selected_functions.wilcox1)[selected_functions.wilcox1$p.value< 0.05]


#Lam and T10
selected_functions2 <- selected_functions[,c(idx_Lam, idx_T10)]
clin.sub_selected_functions2 <-  data.frame(ID = 1:ncol(selected_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_selected_functions2$group)
  return(res$p.value)
}

selected_functions.wilcox2 <-apply(as.matrix(selected_functions2), 1, fx)
selected_functions.wilcox2 <- data.frame(selected_functions.wilcox2 ) 
colnames(selected_functions.wilcox2) <- c("p.value")
write.csv(selected_functions.wilcox2, file = "selected_functions.wilcox0.05_Lam_T10.csv")
rownames(selected_functions.wilcox2)[selected_functions.wilcox2$p.value< 0.01]
rownames(selected_functions.wilcox2)[selected_functions.wilcox2$p.value< 0.05]


#T10 and T4
selected_functions3 <-selected_functions[,c(idx_T4, idx_T10)]
clin.sub_selected_functions3 <-  data.frame(ID = 1:ncol(selected_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_selected_functions3$group)
  return(res$p.value)
}

selected_functions.wilcox3 <-apply(as.matrix(selected_functions3), 1, fx)
selected_functions.wilcox3 <- data.frame(selected_functions.wilcox3 ) 
colnames(selected_functions.wilcox3) <- c("p.value")
write.csv(selected_functions.wilcox3, file = "selected_functions.wilcox0.05_T4_T10.csv")
rownames(selected_functions.wilcox3)[selected_functions.wilcox3$p.value< 0.05]
rownames(selected_functions.wilcox3)[selected_functions.wilcox3$p.value< 0.05]

```


```{r}
selected_functions=read.csv("selected_functions_MAGs.csv",header=T,row.names=NULL,sep=",")
row.names(selected_functions)=selected_functions[,1]
selected_functions=selected_functions[,-1]

selected_functions4=t(selected_functions)
selected_functions5=cbind(selected_functions4, Colvec_expedition)

pairwise.wilcox.test(as.numeric(selected_functions5$`choloylglycine hydrolase`),selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`prephenate dehydratase`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`anthranilate synthase component II `), mOTUs_species5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`indole-3-glycerol phosphate synthase `), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`chorismate synthase`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`pyridoxine 4-dehydrogenase`),selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`dihydrofolate reductase`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`dihydroneopterin aldolase `),selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$lactocepin), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`lactose-specific components`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`N-acetylgalactosamine-specific components`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`galactosamine-specific components`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`galactitol-specific components`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`mannose-specific components `), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`mannitol-specific components `), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$` glucitol/sorbitol-specific components`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`endo-1,4-beta-xylanase`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`arabinan endo-1,5-alpha-L-arabinosidase `), selected_functions3$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`mannan endo-1,4-beta-mannosidase`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`maltose phosphorylase`), selected_functions5$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(as.numeric(selected_functions5$`cellobiose phosphorylase`), selected_functions5$Zones,p.adjust.method = "BH")

```

##C. Selected functions which are differentially abundant in different groups
```{r}
#remove FDR>0.05
selected_functions_new=read.csv("selected_functions_MAGs_new.csv",header=T,row.names=NULL,sep=",")
row.names(selected_functions_new)=selected_functions_new[,1]
selected_functions_new=selected_functions_new[,-1]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

selected_functions_norm <- t(apply(selected_functions_new, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(selected_functions_new)
pheatmap(selected_functions_norm,cluster_rows=FALSE,cluster_cols=FALSE,border_color = NA,annotation_col = my_sample_col)
```

#Figure 5. Phage communities are altered after SCI 
##A. Principal coordinate analyses (PCoA) of Bray-Curtis distances showing that viral community are different between Lam, T10 and T4
```{r}
#Read the abundances of viral populations(>5kb)
pOTUs=read.table("Viral_pops_5KB_95-80_RPKM_only.txt",header=T,row.names=1,sep="\t")
# Transpose the table, phage OTUs
pOTUs1=t(pOTUs)
##a, Principal component analysis showing that phage community abundances are different between healthy and spinal cord injury animals.
# transformate data first
d.clr = log2(pOTUs1+1)
# Calculating Bray-Curtis dissimilarity
Bray=vegdist(d.clr,"bray")
# anosim
ano=anosim(Bray, permutations = 999, grouping=Colvec_expedition$Zones)
summary(ano)
#ANOSIM statistic R: 0.3769 , Significance p: 0.001 

# Applying a non-constrained ordination (PCoA) on the dissimilarity matrix 
pOTUs_bray=capscale(d.clr~-1)
# Extracting the percentage explained by the first two dimensions and automatically adding them to the axes titles
pOTUs_bray_eig = eigenvals(pOTUs_bray)
percentage_variance_explained <- pOTUs_bray_eig / sum(pOTUs_bray_eig)
sum_percentage_variance_explained <- cumsum(pOTUs_bray_eig / sum(pOTUs_bray_eig))
xlabel= as.numeric(format(round((percentage_variance_explained[1]*100), 2), nsmall = 2))
xlabel= sprintf("%.2f %%", xlabel)
xlabel= paste ("PCo1 (", xlabel, "of variation explained", ")")
ylabel= as.numeric(format(round((percentage_variance_explained[2]*100), 2), nsmall = 2))
ylabel= sprintf("%.2f %%", ylabel)
ylabel= paste ("PCo2 (", ylabel,"of variation explained", ")")
# Applying hierarchical clustering on the dissimilarity matrix for plotting on top of the ordination analysis
#H_CLustering=hclust(vegdist(d.clr,"bray"))
# Adding the axes, grid, and other aestethics
#pdf("plots.pdf",width=7, height=7)
plot(pOTUs_bray, type="n",xlab="", ylab="",cex.axis=1, tck = -0.01, mgp = c(3, 1, 0),xlim=c(-5,5),
     xaxp  = c(-1, 1, 2), panel.first=grid(col = "white",lty=0))
axis(1, at=seq(-5, 5, by=1))
axis(2,at=seq(-4,5,by=1))
title(ylab=ylabel, line=2, cex.lab=1)
title(xlab=xlabel, line=2, cex.lab=1)
abline(h=0, v=0, col = "white", lty = 1, lwd = 1.5)
abline(h=-10:10, v=-10:10, col = "lightgray", lty = "dotted", lwd = 0.5)
par(lty=2)
# Adding the hierarchical clustering
colvec_all <- c("#999999", "#999999","#999999","#999999","#999999","#E69F00","#E69F00","#E69F00","#E69F00","#E69F00", "#56B4E9","#56B4E9","#56B4E9","#56B4E9","#56B4E9")

points(pOTUs_bray,cex= 2.5,pch=21, col="black", bg= colvec_all, lwd = 1)
# Adding a legend, bty  The allowed values are "o" (the default) and "n".
legend(-5.5,5.5, pt.cex= 2.5 , pt.lwd = 1,c("Lam", "T10", "T4"),bty = "n" ,pch = 21,col="black",pt.bg = c("#999999", "#E69F00", "#56B4E9"), cex = 1.5,box.lty=2)
ordihull(pOTUs_bray, Colvec_expedition$Zones, display = "sites",label = T)
# ANOSIM, p=0.002. Each data point indicates an individual mouse sample (n = 15).
```

##B.  Whthin-group Bray_Curtis dissimilarity of the viral communities in the Lam, T4, and T10
```{r}
Bray_pOTUs=read.csv("bray_pOTUs.csv",header=T,sep=",")

#Wilcoxon rank sum test corrected by BH
pairwise.wilcox.test(Bray_pOTUs$Bray_curtis, Bray_pOTUs$Zones,p.adjust.method = "BH")

boxplot_Bray_pOTUs=ggplot(Bray_pOTUs,  aes(x=Zones, y=Bray_curtis)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title=NULL ,x=NULL,y = "Bray–Curtis dissimilarity")+theme_classic()
boxplot_Bray_pOTUs1=boxplot_Bray_pOTUs+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
boxplot_Bray_pOTUs1+scale_y_continuous(expand = c(0,0), limits = c(0.1,0.75))+geom_segment(y=0.65,yend=0.65,x=1,xend=3)+geom_text(x=2,y=0.65,label="***",size=10)+geom_segment(y=0.6,yend=0.6,x=1,xend=2)+geom_text(x=1.5,y=0.6,label="**",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```
Note: A higher score suggests higher dissimilarity of different samples in the same group. 

##Between-group Bray_Curtis dissimilarity 
##Bray_Curtis dissimalarity between the control group and SCI groups
```{r}
Bray_pOTUs_SCI=read.csv("Bray_pOTUs_SCI.csv",header=T,sep=",")
wilcox.test(Bray_pOTUs_SCI$Bray_curtis~Bray_pOTUs_SCI$Zones)
pairwise.wilcox.test(Bray_pOTUs_SCI$Bray_curtis,Bray_pOTUs_SCI$Zones,p.adjust.method = "fdr")
#p-value = 0.13
boxplot_Bray_pOTUs_SCI=ggplot(Bray_pOTUs_SCI,  aes(x=Zones, y=Bray_curtis)) + geom_boxplot(fill=c( "#E69F00", "#56B4E9"))+geom_jitter(width = 0.06)+labs(title=NULL,x=NULL,y =c( "Bray_curtis dissimilarity"))+theme_classic()+geom_point(shape=16,fill="black",size=2)

boxplot_Bray_pOTUs_SCI1=boxplot_Bray_pOTUs_SCI+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20,face = "bold"))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
#Including boxplots
p=boxplot_Bray_pOTUs_SCI1+scale_y_continuous(expand = c(0,0), limits=c(0.3,0.7))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
p

```


```{r}
 # Calculating Shannon's "H"
pOTUs_Shannon=diversity(pOTUs1,"shannon")
pOTUs_Shannon=as.data.frame(pOTUs_Shannon)
colnames(pOTUs_Shannon)=c("Shannon")
pOTUs_Shannon1=cbind(pOTUs_Shannon,Colvec_expedition)
#Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(pOTUs_Shannon1$Shannon, pOTUs_Shannon1$Zones,p.adjust.method = "BH")
```

##C.Shannon’s H of the viral communities between the Lam, T4, and T10.
```{r}
#Draw Boxplot
Boxplot_pOTUs_Shannon=ggplot(pOTUs_Shannon1, aes(x=Zones, y=Shannon)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Shannon's H")+theme_classic()
Boxplot_pOTUs_Shannon1=Boxplot_pOTUs_Shannon+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Boxplot_pOTUs_Shannon1+scale_y_continuous(expand = c(0,0), limits = c(5,7.2))+geom_segment(y=7,yend=7,x=1,xend=3)+geom_text(x=2,y=7,label="*",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

```


```{r}
#Differential analysis for viral populations(>5kb) in different treatments 
log2.pOTUs <- log2(pOTUs+1)
#Differential expression for Lam and T10
set.seed(614)                 
log2.pOTUs1 <- log2.pOTUs[,c(idx_Lam, idx_T10)]
clin.sub1 <-  data.frame(ID = 1:ncol(log2.pOTUs1), group = rep(c("Lam","T10"), each=5))

design1<-model.matrix(~ 0 + group, data = clin.sub1)
#head(design1)
contrast.matrix1 <- makeContrasts(groupLam - groupT10, levels=design1)
dge1 <- DGEList(counts=log2.pOTUs1, group = clin.sub1$group) 
v1<-voom(dge1, design1, span = 0.5, plot = TRUE)
##Fitting the linear model
fit1 <- lmFit(v1, design1)
fitc1 <- eBayes(contrasts.fit(fit1, contrast.matrix1))

pOTUs.tt.group1<- topTable(fitc1, coef=1, number = nrow(v1), adjust.method = "BH")
sum(pOTUs.tt.group1$adj.P.Val < 0.05) 

pOTUs.tt.group.voom1 <- rownames(pOTUs.tt.group1)[pOTUs.tt.group1$adj.P.Val< 0.05]
length(pOTUs.tt.group.voom1)
##21

##Differential expression for Lam and T4
idx_Lam <- sample(which(clin.data$group=="Lam"), 5, replace=FALSE) 
idx_T4 <- sample(which(clin.data$group=="T4"), 5, replace=FALSE)
log2.pOTUs2 <- log2.pOTUs[,c(idx_Lam, idx_T4)]
clin.sub2 <-  data.frame(ID = 1:ncol(log2.pOTUs2), group = rep(c("Lam","T4"), each=5))

design2<-model.matrix(~ 0 + group, data = clin.sub2)
#head(design2)
contrast.matrix2 <- makeContrasts(groupLam - groupT4, levels=design2)
dge2 <- DGEList(counts=log2.pOTUs2, group = clin.sub2$group) 
v2<-voom(dge2, design2, span = 0.5, plot = TRUE)
##Fitting the linear model
#double check them whether we can use them in the abundance dataset

fit2 <- lmFit(v2, design2)
fitc2 <- eBayes(contrasts.fit(fit2, contrast.matrix2))

pOTUs.tt.group2<- topTable(fitc2, coef=1, number = nrow(v2), adjust.method = "BH")
sum(pOTUs.tt.group2$adj.P.Val < 0.05)

pOTUs.tt.group.voom2 <- rownames(pOTUs.tt.group2)[pOTUs.tt.group2$adj.P.Val< 0.05]
length(pOTUs.tt.group.voom2)
##57


##Differential expression for T4 and T10
idx_T10 <- sample(which(clin.data$group=="T10"), 5, replace=FALSE) 
idx_T4 <- sample(which(clin.data$group=="T4"), 5, replace=FALSE)
log2.pOTUs3 <- log2.pOTUs[,c(idx_T10, idx_T4)]
clin.sub3 <-  data.frame(ID = 1:ncol(log2.pOTUs3), group = rep(c("T10","T4"), each=5))

design3<-model.matrix(~ 0 + group, data = clin.sub3)
#head(design3)
contrast.matrix3 <- makeContrasts(groupT10 - groupT4, levels=design3)
dge3 <- DGEList(counts=log2.pOTUs3, group = clin.sub3$group) 
v3<-voom(dge3, design3, span = 0.5, plot = TRUE)
##Fitting the linear model
fit3 <- lmFit(v3, design3)
fitc3 <- eBayes(contrasts.fit(fit3, contrast.matrix3))

pOTUs.tt.group3<- topTable(fitc3, coef=1, number = nrow(v3), adjust.method = "BH")
sum(pOTUs.tt.group3$adj.P.Val < 0.05)

pOTUs.tt.group.voom3 <- rownames(pOTUs.tt.group3)[pOTUs.tt.group3$adj.P.Val< 0.05]
length(pOTUs.tt.group.voom3)
##0
```

##D. Volcano plots of t-tests corrected by the Benjamini and Hochberg method for changes to viral populations abundances after spinal cord injury (SCI). 
```{r}
## Volcano plot of Lam and T10 
idx.sig1 <- which(pOTUs.tt.group1$adj.P.Val < 0.05& pOTUs.tt.group1$logFC>0 ) 
idx.sig1_N<-which(pOTUs.tt.group1$adj.P.Val < 0.05& pOTUs.tt.group1$logFC<0 ) 
plot(pOTUs.tt.group1$logFC, -log10(pOTUs.tt.group1$adj.P.Val), pch = 20,xlab = "Difference", ylab = "-log10(P-value)",cex = 1.5,ylim=c(0,3),xlim=c(-2.5,2.5))
points(pOTUs.tt.group1$logFC[idx.sig1], -log10(pOTUs.tt.group1$adj.P.Val)[idx.sig1], pch = 20, col = 2,cex = 1.5)
points(pOTUs.tt.group1$logFC[idx.sig1_N], -log10(pOTUs.tt.group1$adj.P.Val)[idx.sig1_N], pch = 20, col = 4,cex = 1.5)

## Volcano plot of Lam and T4 
idx.sig2 <- which(pOTUs.tt.group2$adj.P.Val < 0.05 & pOTUs.tt.group2$logFC>0 )
idx.sig2_N <- which(pOTUs.tt.group2$adj.P.Val < 0.05 & pOTUs.tt.group2$logFC<0 )
plot(pOTUs.tt.group2$logFC, -log10(pOTUs.tt.group2$adj.P.Val), pch = 20,xlab = "Difference", ylab = "-log10(P-value)",cex = 1.5,ylim=c(0,3),xlim=c(-2.5,2.5))
points(pOTUs.tt.group2$logFC[idx.sig2], -log10(pOTUs.tt.group2$adj.P.Val)[idx.sig2], pch = 20, col = 2,cex = 1.5)
points(pOTUs.tt.group2$logFC[idx.sig2_N], -log10(pOTUs.tt.group2$adj.P.Val)[idx.sig2_N], pch = 20, col = 4,cex = 1.5)

## Volcano plot of T10 and T4 
idx.sig3 <- which(pOTUs.tt.group3$adj.P.Val < 0.05) 
plot(pOTUs.tt.group3$logFC, -log10(pOTUs.tt.group3$adj.P.Val), pch = 20,xlab = "Difference", ylab = "-log10(P-value)",cex = 1.5,ylim=c(0,3),xlim=c(-2.5,2.5))
points(pOTUs.tt.group3$logFC[idx.sig3], -log10(pOTUs.tt.group3$adj.P.Val)[idx.sig3], pch = 20, col = 2)
```


```{r}
predicted_pOTUs_genus=read.csv("predicted_hosts_genus_RPKM.csv",header=T,sep=",")
rownames(predicted_pOTUs_genus)=predicted_pOTUs_genus[,1]
predicted_pOTUs_genus1=predicted_pOTUs_genus[,-1]
#Differential analysis to look for significant predicted_pOTUs_genus
log2.pOTUs_genus <- log2(predicted_pOTUs_genus1+1)

#Lam and T4
log2.pOTUs_genus1 <- log2.pOTUs_genus[,c(idx_Lam, idx_T4)]
clin.sub_pOTUs_genus1 <-  data.frame(ID = 1:ncol(log2.pOTUs_genus1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_pOTUs_genus1$group)
  return(res$p.value)
}
pOTUs_genus.wilcox1 <-apply(as.matrix(log2.pOTUs_genus1), 1, fx)
pOTUs_genus.wilcox1 <- data.frame(pOTUs_genus.wilcox1 ) 
colnames(pOTUs_genus.wilcox1) <- c("p.value")
write.csv(pOTUs_genus.wilcox1,file="pOTUs_genus.wilcox_Lam_T4.csv")
rownames(pOTUs_genus.wilcox1)[pOTUs_genus.wilcox1$p.value< 0.01]
#**, "CAG-1031" ; "Lactobacillus";"Turicibacter" 
rownames(pOTUs_genus.wilcox1)[pOTUs_genus.wilcox1$p.value< 0.05]
#*, Lactococcus;Acutalibacter; Weissella; UBA7182; CAG-56; ASF356;Clostridium_AJ; 

##Lam and T10
log2.pOTUs_genus2 <- log2.pOTUs_genus[,c(idx_Lam, idx_T10)]
clin.sub_pOTUs_genus2 <-  data.frame(ID = 1:ncol(log2.pOTUs_genus2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_pOTUs_genus2$group)
  return(res$p.value)
}
pOTUs_genus.wilcox2 <-apply(as.matrix(log2.pOTUs_genus2), 1, fx)
pOTUs_genus.wilcox2 <- data.frame(pOTUs_genus.wilcox2) 
colnames(pOTUs_genus.wilcox2) <- c("p.value")
write.csv(pOTUs_genus.wilcox2,file="pOTUs_genus.wilcox_Lam_T10.csv")
rownames(pOTUs_genus.wilcox2)[pOTUs_genus.wilcox2$p.value< 0.01]
#**, Eubacterium_R; UBA7160;ASF356
rownames(pOTUs_genus.wilcox2)[pOTUs_genus.wilcox2$p.value< 0.05]
#*, Lactobacillus; UBA3282; Weissella

####T4 and T10
log2.pOTUs_genus3 <- log2.pOTUs_genus[,c(idx_T4, idx_T10)]
clin.sub_pOTUs_genus3 <-  data.frame(ID = 1:ncol(log2.pOTUs_genus3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_pOTUs_genus3$group)
  return(res$p.value)
}
pOTUs_genus.wilcox3 <-apply(as.matrix(log2.pOTUs_genus3), 1, fx)
pOTUs_genus.wilcox3 <- data.frame(pOTUs_genus.wilcox3) 
colnames(pOTUs_genus.wilcox3) <- c("p.value")
write.csv(pOTUs_genus.wilcox3,file="pOTUs_genus.wilcox_T4_T10.csv")
rownames(pOTUs_genus.wilcox3)[pOTUs_genus.wilcox3$p.value< 0.01]
#**,ASF356
rownames(pOTUs_genus.wilcox3)[pOTUs_genus.wilcox3$p.value< 0.05]
#*, Lactococcus, Weissella, Clostridium
```

```{r}
predicted_pOTUs_genus3=t(predicted_pOTUs_genus1)
predicted_pOTUs_genus3=as.data.frame(predicted_pOTUs_genus3)
predicted_pOTUs_genus4=cbind(predicted_pOTUs_genus3,Colvec_expedition)

#Wilcoxon Rank Sum test with BH correction
#p-value<0.05, FDR<0.05
pairwise.wilcox.test(predicted_pOTUs_genus4$`g__CAG-1031`, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Lactobacillus, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Turicibacter, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Weissella, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Lactococcus, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")

#class Clostridia, p-value<0.05, FDR<0.1
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Acutalibacter, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__UBA7182, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__ASF356, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Clostridium_AJ, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Eubacterium_R, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__UBA7160, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$g__UBA3282, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
pairwise.wilcox.test(predicted_pOTUs_genus4$`g__CAG-56`, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")

#others
pairwise.wilcox.test(predicted_pOTUs_genus4$g__Clostridium, predicted_pOTUs_genus4$Zones,p.adjust.method = "BH")
#FDR>0.1

```

#Figure 6. Viral host-prediction reveals that phage abundances vary with their hosts 
## A.phages with FDR<0.05
```{r}
#Only keep phages FDR<0.05
predicted_pOTUs_FDR_0.05=data_frame(predicted_pOTUs_genus4$g__Acutalibacter,predicted_pOTUs_genus4$g__ASF356,predicted_pOTUs_genus4$g__UBA7160,predicted_pOTUs_genus4$g__UBA3282,predicted_pOTUs_genus4$g__Eubacterium_R,predicted_pOTUs_genus4$`g__CAG-1031`,predicted_pOTUs_genus4$g__Lactobacillus,predicted_pOTUs_genus4$g__Turicibacter,predicted_pOTUs_genus4$g__Weissella,predicted_pOTUs_genus4$g__Lactococcus)

colnames(predicted_pOTUs_FDR_0.05)=c("Acutalibacter","ASF356","UBA7160","UBA3282","Eubacterium_R","CAG-1031","Lactobacillus","Turicibacter","Weissella","Lactococcus")
row.names(predicted_pOTUs_FDR_0.05)=row.names(singleM_genus)
predicted_pOTUs_FDR_0.05_t=t(predicted_pOTUs_FDR_0.05)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

predicted_pOTUs_FDR_0.05_norm <- t(apply(predicted_pOTUs_FDR_0.05_t, 1, cal_z_score))

pheatmap(predicted_pOTUs_FDR_0.05_norm,cluster_cols=FALSE, legend = TRUE,annotation_col = my_sample_col,border_color = NA)


```

##B-D. Phages that infect CAG-1031, Lactobacillus and Turicibacter decreased after spinal cord injury
```{r}
#CAG_1031_phage
CAG_1031_phage=ggplot(predicted_pOTUs_genus4, aes(x=Zones, y=`g__CAG-1031`))+ geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "CAG-1031_phage",x=NULL,y = "RPKM")+theme_classic()

CAG_1031_boxplot=CAG_1031_phage+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
CAG_1031_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,750))+geom_segment(y=700,yend=700,x=1,xend=3)+geom_text(x=2,y=700,label="**",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

#Lactobacillus_phage
Lactobacillus_phage=ggplot(predicted_pOTUs_genus4, aes(x=Zones, y=g__Lactobacillus))+ geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "Lactobacillus_phage",x=NULL,y = "RPKM")+theme_classic()
Lactobacillus_boxplot=Lactobacillus_phage+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Lactobacillus_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,400))+geom_segment(y=350,yend=350,x=1,xend=2)+geom_segment(y=380,yend=380,x=1,xend=3)+geom_text(x=1.5,y=350,label="*",size=10)+geom_text(x=2,y=380,label="**",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

#Turicibacter
Turicibacter_phage=ggplot(predicted_pOTUs_genus4, aes(x=Zones, y=g__Turicibacter))+ geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "Turicibacter_phage",x=NULL,y = "RPKM")+theme_classic()

Turicibacter_boxplot=Turicibacter_phage+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Turicibacter_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,22))+geom_segment(y=20,yend=20,x=1,xend=3)+geom_text(x=2,y=20,label="**",size=10)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
```

##E,F. Phages that infect Weissella and Lactococcus increased after spinal cord injury.
```{r}
 #Weissella_phage
Weissella_phage=ggplot(predicted_pOTUs_genus4, aes(x=Zones, y=g__Weissella))+ geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "Weissella_phage",x=NULL,y = "RPKM")+theme_classic()

Weissella_boxplot=Weissella_phage+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Weissella_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,350))+geom_segment(y=330,yend=330,x=1,xend=3)+geom_text(x=2,y=330,label="*",size=10)+geom_segment(y=300,yend=300,x=2,xend=3)+geom_text(x=2.5,y=300,label="*",size=10)+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

#Lactococcus_phage
Lactococcus_phage=ggplot(predicted_pOTUs_genus4, aes(x=Zones, y=g__Lactococcus))+ geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "Lactococcus_phage",x=NULL,y = "RPKM")+theme_classic()

Lactococcus_boxplot=Lactococcus_phage+theme(plot.title = element_text(size = 24,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Lactococcus_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,25))+geom_segment(y=23,yend=23,x=1,xend=3)+geom_text(x=2,y=23,label="*",size=10)+geom_segment(y=21,yend=21,x=2,xend=3)+geom_text(x=2.5,y=21,label="*",size=10)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
#The central line in each box represents the median value of the data set. 
```

##G. Phages that infect members in class Clostridia were altered after spinal cord injury.
```{r}
predicted_pOTUs_Clostridia=data_frame(predicted_pOTUs_genus4$g__Acutalibacter,predicted_pOTUs_genus4$g__ASF356,predicted_pOTUs_genus4$g__UBA7160,predicted_pOTUs_genus4$g__UBA3282,predicted_pOTUs_genus4$g__Eubacterium_R,singleM_genus$Zones)

colnames(predicted_pOTUs_Clostridia)=c("Acutalibacter","ASF356","UBA7160","UBA3282","Eubacterium_R","Zones")
row.names(predicted_pOTUs_Clostridia)=row.names(singleM_genus)
predicted_pOTUs_Clostridia1=gather(predicted_pOTUs_Clostridia,key="Clostridia",value = "abundances",-Zones)

predicted_pOTUs_Clostridia1$Clostridia=factor(predicted_pOTUs_Clostridia1$Clostridia,levels=c("ASF356","UBA7160","UBA3282","Acutalibacter","Eubacterium_R"))

pOTUs_Clostridia_boxplot=ggplot(predicted_pOTUs_Clostridia1, aes(x=Clostridia, y=abundances,fill=factor(Zones)))+theme_classic()+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+geom_boxplot(position=position_dodge(0.8))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8),binpositions = "bygroup", dotsize=0.4)+labs(title =NULL,x=NULL,y = "Class Clostridia \n relative abundances(%)")+theme_classic()

pOTUs_Clostridia_boxplot=pOTUs_Clostridia_boxplot+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 10,angle=45,hjust = 1))+theme(axis.text.y = element_text(color="black",size = 10),axis.title.y=element_text(size=12))

pOTUs_Clostridia_boxplot+scale_y_continuous(expand = c(0,0),limits = c(0,100))+geom_text(x=1.275,y=5,label="*",size=8)+geom_text(x=2,y=2,label="**",size=8)+geom_text(x=2.275,y=3,label="*",size=8)+geom_text(x=3,y=15,label="**",size=8)+geom_text(x=3.275,y=10,label="*",size=8)+geom_text(x=4.275,y=19,label="*",size=8)+geom_text(x=5.275,y=70,label="*",size=8)+geom_text(x=6,y=15,label="*",size=8)+geom_text(x=7.275,y=45,label="*",size=8)+geom_text(x=8,y=85,label="**",size=8)+ theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

```



#Supplemental figures

#Figure S2. Differential abundance analysis of bacteria across three treatment groups

##A. Venn diagram of the number of shared and unshared bacterial clustered OTUs in different groups. 
```{r}
library(VennDiagram)
L2_Venn=read.csv("L2_Venn_Diagrams.csv",header=T,row.names=NULL,sep=",")
rownames(L2_Venn)=L2_Venn[,1]
L2_Venn=L2_Venn[,-1]
grid.newpage()
draw.triple.venn(asp=1, area1 = nrow(subset(L2_Venn, Lam == 1)), area2 = nrow(subset(L2_Venn, T10 == 1)), area3 = nrow(subset(L2_Venn, T4 == 1)), n12 = nrow(subset(L2_Venn, Lam == 1 & 
T10== 1)), n23 = nrow(subset(L2_Venn, T10 == 1 & T4 == 1)), n13 = nrow(subset(L2_Venn, Lam == 1 & T4 == 1)), n123 = nrow(subset(L2_Venn, Lam == 1 & T10 == 1 & T4 == 1)), 
category = c("Lam", "T10", "T4"), lty = "blank",fill = c("#999999", "#E69F00", "#56B4E9"))

```

##B. Shannon’s H of the microbial communities between the Lam, T4, and T10.
```{r}
# Calculating Shannon's "H"
L2_Shannon=diversity(L2,"shannon")
L2_Shannon=as.data.frame(L2_Shannon)
colnames(L2_Shannon)=c("Shannon")
L2_Shannon1=cbind(L2_Shannon,Colvec_expedition)
#Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(L2_Shannon1$Shannon, L2_Shannon1$Zones,p.adjust.method = "BH")

#calculate richness
L2_richness <- specnumber(L2)
L2_richness=as.data.frame(L2_richness)
colnames(L2_richness)=c("Richness")
L2_richness1=cbind(L2_richness,Colvec_expedition)
#Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(L2_richness1$Richness, L2_richness1$Zones,p.adjust.method = "BH")

#calculate eveness 
L2_evenness <- L2_Shannon/log(L2_richness)
L2_evenness=as.data.frame(L2_evenness)
colnames(L2_evenness)=c("Evenness")
L2_evenness1=cbind(L2_evenness,Colvec_expedition)
#Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(L2_evenness1$Evenness, L2_evenness1$Zones,p.adjust.method = "BH")

#Combine theree characters
L2_diversity=data.frame(L2_Shannon$Shannon,L2_richness$Richness,L2_evenness$Evenness)

#Draw Boxplot
Boxplot_L2_Shannon=ggplot(L2_Shannon1, aes(x=Zones, y=Shannon)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Shannon's H")+theme_classic()
Boxplot_L2_Shannon1=Boxplot_L2_Shannon+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Boxplot_L2_Shannon1+scale_y_continuous(expand = c(0,0), limits = c(2.8,4.5))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

```


##C.Differential abundances of bacterial genera 
```{r  fig.height = 12, fig.width = 5, fig.align = "center"}
#B.Differential abundances of bacterial genera (p<0.05 by Wilcoxon signed-rank test) in either two groups are indicated in red. Each row representing a unique bacterial genus was Z-score normalized. Bacterial genera on the y-axis are clustered using Euclidean distances.  
#use pheatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

genus_singleM_norm <- t(apply(genus_singleM, 1, cal_z_score))

pheatmap(genus_singleM_norm,cluster_cols=FALSE, legend = TRUE,annotation_col = my_sample_col)

```

##D rare genera relative abundance(%) FDR<0.05
```{r}

singleM_genus_rare=data_frame(singleM_genus$Clostridium_M,singleM_genus$Eubacterium_I,singleM_genus$Slackia,singleM_genus$`CAG-791`, singleM_genus$Flavonifractor,singleM_genus$Zones)
rownames(singleM_genus_rare)=rownames(singleM_genus)

colnames(singleM_genus_rare)=c("Clostridium_M","Eubacterium_I","Slackia","CAG-791","Flavonifractor","Zones")

singleM_genus_rare1=gather(singleM_genus_rare,key="Rare",value = "abundances",-Zones)

singleM_genus_rare1$Rare=factor(singleM_genus_rare1$Rare,levels=c("Clostridium_M","Eubacterium_I","Flavonifractor","Slackia","CAG-791"))

genus_rare_boxplot=ggplot(singleM_genus_rare1, aes(x=Rare, y=abundances,fill=factor(Zones)))+theme_classic()+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+geom_boxplot(position=position_dodge(0.8))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8),binpositions = "bygroup", dotsize=0.4)+labs(title =NULL,x=NULL,y = "Class Clostridia \n relative abundances(%)")+theme_classic()

genus_rare_boxplot=genus_rare_boxplot+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 10,angle=45,hjust = 1))+theme(axis.text.y = element_text(color="black",size = 10),axis.title.y=element_text(size=12))

genus_rare_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,0.55))+geom_segment(y=0.5,yend=0.5,x=0.75,xend=1)+geom_text(x=0.875,y=0.5,label="*",size=6)+geom_segment(y=0.29,yend=0.29,x=1,xend=1.25)+geom_text(x=1.125,y=0.29,label="*",size=6)+geom_segment(y=0.31,yend=0.31,x=2,xend=2.25)+geom_text(x=2.125,y=0.31,label="**",size=6)+geom_segment(y=0.23,yend=0.23,x=2.75,xend=3)+geom_text(x=2.875,y=0.23,label="*",size=6)+geom_segment(y=0.11,yend=0.11,x=3.75,xend=4)+geom_text(x=3.875,y=0.11,label="*",size=6)+geom_segment(y=0.11,yend=0.11,x=5,xend=5.25)+geom_text(x=5.125,y=0.11,label="**",size=6)+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

```


#Figure S3. Species-level differential abundance analysis of bacteria across three treatment groups.
```{r , fig.height = 12, fig.width = 5, fig.align = "center"}
#Differential abundances of bacterial species (p<0.05 by Wilcoxon signed-rank test) in either two groups are indicated in red. Each row representing a unique bacterial genus was Z-score normalized.

mOTUs_species=read.table("MAGs_species.txt",header=T,row.names=1,sep="\t")
mOTUs_species$Maxbins=as.factor(rownames((mOTUs_species)))
mOTUs_species1=merge(mOTUs_species,mOTUs,by="Maxbins")
mOTUs_species2=mOTUs_species1[,-1]
row.names(mOTUs_species2)=mOTUs_species2[,1]
mOTUs_species3=mOTUs_species2[,-1]

mOTUs_species_norm <- t(apply(mOTUs_species3, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(mOTUs_species3)

pheatmap(mOTUs_species_norm,cluster_cols=FALSE, annotation_col = my_sample_col,legend = TRUE)

```


```{r}
# Calculating Shannon's "H"
unique_PCs_Shannon=diversity(t_unique_PCs,"shannon")
unique_PCs_Shannon=as.data.frame(unique_PCs_Shannon)
colnames(unique_PCs_Shannon)=c("Shannon")
unique_PCs_Shannon1=cbind(unique_PCs_Shannon,Colvec_expedition)
#Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(unique_PCs_Shannon1$Shannon, unique_PCs_Shannon1$Zones,p.adjust.method = "BH")
```
#Figure S4. Predicted metabolic pathways are different between Lam controls and SCI groups. 
##A.Shannon’s H of the microbial functions between the Lam, T4, and T10.
```{r}
#Draw Boxplot
Boxplot_unique_PCs_Shannon=ggplot(unique_PCs_Shannon1, aes(x=Zones, y=Shannon)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = NULL,x=NULL,y = "Shannon's H")+theme_classic()


Boxplot_unique_PCs_Shannon1=Boxplot_unique_PCs_Shannon+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))

Boxplot_unique_PCs_Shannon1+scale_y_continuous(expand = c(0,0), limits = c(9.5,10.5))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

```


```{r}
##carbohydrate metabolism
KEGG_carbohydrate=read.csv("mean_KEGGFUN_wilcox_0.05_carbohydrate_metabolism.csv",header=T,sep=",")
rownames(KEGG_carbohydrate)=KEGG_carbohydrate[,1]
KEGG_carbohydrate1=KEGG_carbohydrate[,-(1:2)]
```

```{r}
#Statistics analysis
#Lam and T4
KEGG_carbohydrate_functions1 <- KEGG_carbohydrate1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_carbohydrate_functions1 <-  data.frame(ID = 1:ncol(KEGG_carbohydrate_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_carbohydrate_functions1$group)
  return(res$p.value)
}

KEGG_carbohydrate_functions.wilcox1 <-apply(as.matrix(KEGG_carbohydrate_functions1), 1, fx)
KEGG_carbohydrate_functions.wilcox1 <- data.frame(KEGG_carbohydrate_functions.wilcox1 ) 
colnames(KEGG_carbohydrate_functions.wilcox1) <- c("p.value")
write.csv(KEGG_carbohydrate_functions.wilcox1, file = "KEGG_carbohydrate_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_carbohydrate_functions2 <- KEGG_carbohydrate1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_carbohydrate_functions2 <-  data.frame(ID = 1:ncol(KEGG_carbohydrate_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_carbohydrate_functions2$group)
  return(res$p.value)
}

KEGG_carbohydrate_functions.wilcox2 <-apply(as.matrix(KEGG_carbohydrate_functions2), 1, fx)
KEGG_carbohydrate_functions.wilcox2 <- data.frame(KEGG_carbohydrate_functions.wilcox2 ) 
colnames(KEGG_carbohydrate_functions.wilcox2) <- c("p.value")
write.csv(KEGG_carbohydrate_functions.wilcox2, file = "KEGG_carbohydrate_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_carbohydrate_functions3 <-KEGG_carbohydrate1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_carbohydrate_functions3 <-  data.frame(ID = 1:ncol(KEGG_carbohydrate_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_carbohydrate_functions3$group)
  return(res$p.value)
}
KEGG_carbohydrate_functions.wilcox3 <-apply(as.matrix(KEGG_carbohydrate_functions3), 1, fx)
KEGG_carbohydrate_functions.wilcox3 <- data.frame(KEGG_carbohydrate_functions.wilcox3 ) 
colnames(KEGG_carbohydrate_functions.wilcox3) <- c("p.value")
write.csv(KEGG_carbohydrate_functions.wilcox3, file = "KEGG_carbohydrate_functions.wilcox0.05_T4_T10.csv")

```

```{r}
#Wilcoxon Rank Sum test with false discovery rates 
KEGG_carbohydrate2=t(KEGG_carbohydrate1)
KEGG_carbohydrate3=cbind(KEGG_carbohydrate2, Colvec_expedition)

pairwise.wilcox.test(KEGG_carbohydrate3$`chitinase [EC:3.2.1.14]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`N-acetylmuramic acid 6-phosphate etherase [EC:4.2.-.-]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`N-acetylneuraminate synthase [EC:2.5.1.56]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`N-acylglucosamine 2-epimerase [EC:5.1.3.8]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`NAD-dependent deacetylase [EC:3.5.1.-]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(KEGG_carbohydrate3$`UDP-N-acetylmuramate dehydrogenase [EC:1.1.1.158]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(KEGG_carbohydrate3$`acetoin dehydrogenase [EC:1.1.1.5]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`succinate-semialdehyde dehydrogenase (NADP+) [EC:1.2.1.16]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`hydroxymethylglutaryl-CoA synthase [EC:2.3.3.10]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`isocitrate dehydrogenase (NAD+) [EC:1.1.1.41]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`citrate lyase subunit beta [EC:4.1.3.6]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`succinyl-CoA synthetase alpha subunit [EC:6.2.1.5]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`L-iditol 2-dehydrogenase [EC:1.1.1.14]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`L-rhamnose isomerase [EC:5.3.1.14]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`mannan endo-1,4-beta-mannosidase [EC:3.2.1.78]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, glucitol/sorbitol-specific IIA component [EC:2.7.1.69]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, mannitol-specific IIA component [EC:2.7.1.69]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`mannose-6-phosphate isomerase [EC:5.3.1.8]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, mannose-specific IIA component [EC:2.7.1.69]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, mannose-specific IIB component [EC:2.7.1.69]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, mannose-specific IIC component`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`galactose-6-phosphate isomerase [EC:5.3.1.26]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, galactitol-specific IIA component [EC:2.7.1.69]`,KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, galactitol-specific IIB component [EC:2.7.1.69]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, galactitol-specific IIC component`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, galactosamine-specific IIB component [EC:2.7.1.69]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, galactosamine-specific IIC component`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, galactosamine-specific IID component`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, lactose-specific IIA component [EC:2.7.1.69]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`PTS system, N-acetylgalactosamine-specific IIA component`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`glucose-1-phosphatase [EC:3.1.3.10]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`phosphoglycerate mutase [EC:5.4.2.1]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(KEGG_carbohydrate3$`pyruvate dehydrogenase E2 component (dihydrolipoamide`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`phosphoenolpyruvate carboxykinase (GTP) [EC:4.1.1.32]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`phosphoglucomutase [EC:5.4.2.2]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`acetaldehyde dehydrogenase / alcohol dehydrogenase [EC:1.2.1.10`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(KEGG_carbohydrate3$`S-(hydroxymethyl)glutathione dehydrogenase / alcohol dehydrogenase`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`alcohol dehydrogenase [EC:1.1.1.1]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`phosphoglycolate phosphatase [EC:3.1.3.18]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`formate dehydrogenase, alpha subunit [EC:1.2.1.2]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`formyltetrahydrofolate deformylase [EC:3.5.1.10]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`myo-inositol-1(or 4)-monophosphatase [EC:3.1.3.25]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`rhamnulose-1-phosphate aldolase [EC:4.1.2.19]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`3-hexulose-6-phosphate synthase [EC:4.1.2.43]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`glucose 1-dehydrogenase [EC:1.1.1.47]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`ribokinase [EC:2.7.1.15]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`phosphopentomutase [EC:5.4.2.7]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`methylglyoxal synthase [EC:4.2.3.3]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`lactaldehyde reductase [EC:1.1.1.77]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`acetyl-CoA carboxylase biotin carboxyl carrier protein`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`propionate CoA-transferase [EC:2.8.3.1]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`phosphoenolpyruvate carboxylase [EC:4.1.1.31]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`cellobiose phosphorylase [EC:2.4.1.20]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`maltose phosphorylase [EC:2.4.1.8]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`arabinan endo-1,5-alpha-L-arabinosidase [EC:3.2.1.99]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`endo-1,4-beta-xylanase [EC:3.2.1.8]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`myo-inositol catabolism protein IolS [EC:1.1.1.-]`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`putative family 31 glucosidase`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_carbohydrate3$`sugar fermentation stimulation protein A`, KEGG_carbohydrate3$Zones,p.adjust.method = "fdr")
```

##carbohydrate metabolism
```{r , fig.height = 15, fig.width = 10, fig.align = "center"}
#Use FDR<0.05
KEGG_carbohydrate_new=read.csv("mean_KEGGFUN_FDR_0.05_carbohydrate_metabolism.csv",header=T,sep=",")
rownames(KEGG_carbohydrate_new)=KEGG_carbohydrate_new[,1]
KEGG_carbohydrate_new1=KEGG_carbohydrate_new[,-(1:2)]
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
KEGG_carbohydrate_new1_norm <- t(apply(KEGG_carbohydrate_new1, 1, cal_z_score))

my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_carbohydrate_new1_norm)
#pdf(file = "carbohydrate_metabolism", 
   # width = 10, # The width of the plot in inches
   # height = 15)
pheatmap(KEGG_carbohydrate_new1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE, cluster_rows = T,cluster_cols = F,border_color = NA,cellwidth = 12, cellheight = 15, main="Carbohydrate Metabolism")

```


```{r}
##energy_metabolism
KEGG_energy=read.csv("mean_KEGGFUN_wilcox_0.05_energy_metabolism.csv",header=T,sep=",")
rownames(KEGG_energy)=KEGG_energy[,1]
KEGG_energy1=KEGG_energy[,-(1:2)]

```

```{r}
#Lam and T4
KEGG_energy_functions1 <- KEGG_energy1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_energy_functions1 <-  data.frame(ID = 1:ncol(KEGG_energy_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_energy_functions1$group)
  return(res$p.value)
}

KEGG_energy_functions.wilcox1 <-apply(as.matrix(KEGG_energy_functions1), 1, fx)
KEGG_energy_functions.wilcox1 <- data.frame(KEGG_energy_functions.wilcox1 ) 
colnames(KEGG_energy_functions.wilcox1) <- c("p.value")
write.csv(KEGG_energy_functions.wilcox1, file = "KEGG_energy_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_energy_functions2 <- KEGG_energy1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_energy_functions2 <-  data.frame(ID = 1:ncol(KEGG_energy_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_energy_functions2$group)
  return(res$p.value)
}

KEGG_energy_functions.wilcox2 <-apply(as.matrix(KEGG_energy_functions2), 1, fx)
KEGG_energy_functions.wilcox2 <- data.frame(KEGG_energy_functions.wilcox2 ) 
colnames(KEGG_energy_functions.wilcox2) <- c("p.value")
write.csv(KEGG_energy_functions.wilcox2, file = "KEGG_energy_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_energy_functions3 <-KEGG_energy1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_energy_functions3 <-  data.frame(ID = 1:ncol(KEGG_energy_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_energy_functions3$group)
  return(res$p.value)
}
KEGG_energy_functions.wilcox3 <-apply(as.matrix(KEGG_energy_functions3), 1, fx)
KEGG_energy_functions.wilcox3 <- data.frame(KEGG_energy_functions.wilcox3 ) 
colnames(KEGG_energy_functions.wilcox3) <- c("p.value")
write.csv(KEGG_energy_functions.wilcox3, file = "KEGG_energy_functions.wilcox0.05_T4_T10.csv")
```

```{r}
KEGG_energy2=t(KEGG_energy1)
KEGG_energy3=cbind(KEGG_energy2, Colvec_expedition)

pairwise.wilcox.test(KEGG_energy3$`methylenetetrahydrofolate reductase (NADPH) [EC:1.5.1.20]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`carbonic anhydrase [EC:4.2.1.1]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`Nif-specific regulatory protein`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`glutamate synthase (NADPH/NADH) large chain [EC:1.4.1.13 1.4.1.14]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`aspartate--ammonia ligase [EC:6.3.1.1]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`cystathionine beta-lyase [EC:4.4.1.8]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`NADH dehydrogenase I subunit F [EC:1.6.5.3]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`V-type H+-transporting ATPase subunit A [EC:3.6.3.14]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`V-type H+-transporting ATPase subunit E [EC:3.6.3.14]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`F-type H+-transporting ATPase subunit c [EC:3.6.3.14]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`cystathionine gamma-synthase [EC:2.5.1.48]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`adenylylsulfate kinase [EC:2.7.1.25]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`dihydroorotate dehydrogenase electron transfer subunit`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`electron transfer flavoprotein alpha subunit`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`electron transport complex protein RnfG`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`flavodoxin I`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`hydrogenase-4 component B [EC:1.-.-.-]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`NADPH2:quinone reductase [EC:1.6.5.5]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`phenylacetic acid degradation protein`, KEGG_energy3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_energy3$`Cu2+-exporting ATPase [EC:3.6.3.4]`, KEGG_energy3$Zones,p.adjust.method = "fdr")
```

##energy_metabolism
```{r , fig.height = 15, fig.width = 10, fig.align = "center"}
KEGG_energy_new=read.csv("mean_KEGGFUN_FDR_0.05_energy_metabolism.csv",header=T,sep=",")
rownames(KEGG_energy_new)=KEGG_energy_new[,1]
KEGG_energy_new1=KEGG_energy_new[,-(1:2)]

KEGG_energy_new1_norm <- t(apply(KEGG_energy_new1, 1, cal_z_score))

my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_energy_new1_norm)
#pdf(file = "energy_metabolism", 
    #width = 10, # The width of the plot in inches
    #height = 15)

pheatmap(KEGG_energy_new1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE, cluster_rows = T,cluster_cols = F,border_color = NA,cellwidth = 12, cellheight = 15, main="Energy Metabolism")

```


##amino acid metabolism
```{r}
KEGG_amino_acid=read.csv("mean_KEGGFUN_wilcox_0.05_amino.acids.csv",header=T,sep=",")
rownames(KEGG_amino_acid)=KEGG_amino_acid[,1]
KEGG_amino_acid1=KEGG_amino_acid[,-(1:2)]
```


```{r}
#Lam and T4
KEGG_amino_acid_functions1 <- KEGG_amino_acid1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_amino_acid_functions1 <-  data.frame(ID = 1:ncol(KEGG_amino_acid_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_amino_acid_functions1$group)
  return(res$p.value)
}

KEGG_amino_acid_functions.wilcox1 <-apply(as.matrix(KEGG_amino_acid_functions1), 1, fx)
KEGG_amino_acid_functions.wilcox1 <- data.frame(KEGG_amino_acid_functions.wilcox1 ) 
colnames(KEGG_amino_acid_functions.wilcox1) <- c("p.value")
write.csv(KEGG_amino_acid_functions.wilcox1, file = "KEGG_amino_acid_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_amino_acid_functions2 <- KEGG_amino_acid1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_amino_acid_functions2 <-  data.frame(ID = 1:ncol(KEGG_amino_acid_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_amino_acid_functions2$group)
  return(res$p.value)
}

KEGG_amino_acid_functions.wilcox2 <-apply(as.matrix(KEGG_amino_acid_functions2), 1, fx)
KEGG_amino_acid_functions.wilcox2 <- data.frame(KEGG_amino_acid_functions.wilcox2 ) 
colnames(KEGG_amino_acid_functions.wilcox2) <- c("p.value")
write.csv(KEGG_amino_acid_functions.wilcox2, file = "KEGG_amino_acid_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_amino_acid_functions3 <-KEGG_amino_acid1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_amino_acid_functions3 <-  data.frame(ID = 1:ncol(KEGG_amino_acid_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_amino_acid_functions3$group)
  return(res$p.value)
}
KEGG_amino_acid_functions.wilcox3 <-apply(as.matrix(KEGG_amino_acid_functions3), 1, fx)
KEGG_amino_acid_functions.wilcox3 <- data.frame(KEGG_amino_acid_functions.wilcox3 ) 
colnames(KEGG_amino_acid_functions.wilcox3) <- c("p.value")
write.csv(KEGG_amino_acid_functions.wilcox3, file = "KEGG_amino_acid_functions.wilcox0.05_T4_T10.csv")
```

```{r}
#Wilcoxon Rank Sum test with false discovery rates 

KEGG_amino_acid2=t(KEGG_amino_acid1)
KEGG_amino_acid3=cbind(KEGG_amino_acid2, Colvec_expedition)

pairwise.wilcox.test(KEGG_amino_acid3$`glutamate 5-kinase [EC:2.7.2.11]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`glutamate N-acetyltransferase / amino-acid N-acetyltransferase`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`arginine deiminase [EC:3.5.3.6]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`agmatinase [EC:3.5.3.11]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`amidase [EC:3.5.1.4]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`ornithine decarboxylase [EC:4.1.1.17]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`S-adenosylmethionine decarboxylase [EC:4.1.1.50]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`D-serine dehydratase [EC:4.3.1.18]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`threonine synthase [EC:4.2.3.1]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`cyclase HisF [EC:4.1.3.-]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`phosphoribosyl-ATP pyrophosphohydrolase / phosphoribosyl-AMP`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`imidazoleglycerol-phosphate dehydratase / histidinol-phosphatase`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`phosphoribosyl-ATP pyrophosphohydrolase [EC:3.6.1.31]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`phosphoribosyl-AMP cyclohydrolase [EC:3.5.4.19]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`2-aminoadipate transaminase [EC:2.6.1.-]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`D-alanine transaminase [EC:2.6.1.21]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")

pairwise.wilcox.test(KEGG_amino_acid3$`chorismate synthase [EC:4.2.3.5]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`indole-3-glycerol phosphate synthase [EC:4.1.1.48]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`anthranilate synthase component II [EC:4.1.3.27]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`prephenate dehydratase [EC:4.2.1.51]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`D-citramalate synthase [EC:2.3.1.182]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`3-hydroxyisobutyrate dehydrogenase [EC:1.1.1.31]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`leucyl aminopeptidase [EC:3.4.11.1]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`phosphinothricin acetyltransferase [EC:2.3.1.183]`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_amino_acid3$`glyoxylase I family protein`, KEGG_amino_acid3$Zones,p.adjust.method = "fdr")
```

##amino_acid_metabolism
```{r, fig.height = 15, fig.width = 10, fig.align = "center"}
KEGG_amino_acid_new=read.csv("mean_KEGGFUN_FDR_0.05_amino.acids.csv",header=T,sep=",")
rownames(KEGG_amino_acid_new)=KEGG_amino_acid_new[,1]
KEGG_amino_acid_new1=KEGG_amino_acid_new[,-(1:2)]

KEGG_amino_acid_new1_norm <- t(apply(KEGG_amino_acid_new1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_amino_acid_new1_norm)
#pdf(file="Amino Acid Metabolism", 
         #width = 10, # The width of the plot in inches
   # height = 15)
pheatmap(KEGG_amino_acid_new1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE, cluster_rows = T,cluster_cols = F,border_color = NA,cellwidth = 12, cellheight = 15,main="Amino Acid Metabolism")

```


```{r, fig.height = 15, fig.width = 10, fig.align = "center"}
##lipid metabolism
KEGG_lipid=read.csv("mean_KEGGFUN_wilcox_0.05_Lipid.Metabolism.csv",header=T,sep=",")
rownames(KEGG_lipid)=KEGG_lipid[,1]
KEGG_lipid1=KEGG_lipid[,-(1:2)]
```

```{r}
#Lam and T4
KEGG_lipid_functions1 <- KEGG_lipid1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_lipid_functions1 <-  data.frame(ID = 1:ncol(KEGG_lipid_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_lipid_functions1$group)
  return(res$p.value)
}

KEGG_lipid_functions.wilcox1 <-apply(as.matrix(KEGG_lipid_functions1), 1, fx)
KEGG_lipid_functions.wilcox1 <- data.frame(KEGG_lipid_functions.wilcox1 ) 
colnames(KEGG_lipid_functions.wilcox1) <- c("p.value")
write.csv(KEGG_lipid_functions.wilcox1, file = "KEGG_lipid_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_lipid_functions2 <- KEGG_lipid1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_lipid_functions2 <-  data.frame(ID = 1:ncol(KEGG_lipid_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_lipid_functions2$group)
  return(res$p.value)
}

KEGG_lipid_functions.wilcox2 <-apply(as.matrix(KEGG_lipid_functions2), 1, fx)
KEGG_lipid_functions.wilcox2 <- data.frame(KEGG_lipid_functions.wilcox2 ) 
colnames(KEGG_lipid_functions.wilcox2) <- c("p.value")
write.csv(KEGG_lipid_functions.wilcox2, file = "KEGG_lipid_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_lipid_functions3 <-KEGG_lipid1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_lipid_functions3 <-  data.frame(ID = 1:ncol(KEGG_lipid_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_lipid_functions3$group)
  return(res$p.value)
}
KEGG_lipid_functions.wilcox3 <-apply(as.matrix(KEGG_lipid_functions3), 1, fx)
KEGG_lipid_functions.wilcox3 <- data.frame(KEGG_lipid_functions.wilcox3 ) 
colnames(KEGG_lipid_functions.wilcox3) <- c("p.value")
write.csv(KEGG_lipid_functions.wilcox3, file = "KEGG_lipid_functions.wilcox0.05_T4_T10.csv")
```

```{r}
#Wilcoxon Rank Sum test with "BH" (also known as "fdr") 
KEGG_lipid2=t(KEGG_lipid1)
KEGG_lipid3=cbind(KEGG_lipid2, Colvec_expedition)
pairwise.wilcox.test(KEGG_lipid3$`gamma-glutamyltranspeptidase [EC:2.3.2.2]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`acetyl/propionyl carboxylase subunit alpha`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`3R-hydroxymyristoyl ACP dehydrase [EC:4.2.1.-]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`enoyl-[acyl-carrier protein] reductase I [EC:1.3.1.9]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`rubredoxin-NAD+ reductase [EC:1.18.1.1]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`1,3-propanediol dehydrogenase [EC:1.1.1.202]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`dihydroxyacetone kinase, N-terminal domain [EC:2.7.1.-]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`fatty-acyl-CoA synthase [EC:6.2.1.-]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`choloylglycine hydrolase [EC:3.5.1.24]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_lipid3$`arylsulfatase [EC:3.1.6.1]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")

```

##lipid metabolism
```{r , fig.height = 15, fig.width = 10, fig.align = "center"}
##lipid metabolism
KEGG_lipid_new=read.csv("mean_KEGGFUN_FDR_0.05_Lipid.Metabolism.csv",header=T,sep=",")
rownames(KEGG_lipid_new)=KEGG_lipid_new[,1]
KEGG_lipid_new1=KEGG_lipid_new[,-(1:2)]

KEGG_lipid_new1_norm <- t(apply(KEGG_lipid_new1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_lipid_new1_norm)
#pdf(file="Lipid Metabolism",height = 15,width = 10)
pheatmap(KEGG_lipid_new1_norm,show_colnames = F, show_rownames = T,legend = TRUE,border_color = NA,cluster_cols = F,cellwidth = 12, cellheight = 15,main="Lipid Metabolism")
```


```{r}
##enzyme metabolism
KEGG_enzyme=read.csv("mean_KEGGFUN_wilcox_0.05_enzyme.families.csv",header=T,sep=",")
rownames(KEGG_enzyme)=KEGG_enzyme[,1]
KEGG_enzyme1=KEGG_enzyme[,-(1:2)]
```

```{r}
#Lam and T4
KEGG_enzyme_functions1 <- KEGG_enzyme1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_enzyme_functions1 <-  data.frame(ID = 1:ncol(KEGG_enzyme_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_enzyme_functions1$group)
  return(res$p.value)
}

KEGG_enzyme_functions.wilcox1 <-apply(as.matrix(KEGG_enzyme_functions1), 1, fx)
KEGG_enzyme_functions.wilcox1 <- data.frame(KEGG_enzyme_functions.wilcox1 ) 
colnames(KEGG_enzyme_functions.wilcox1) <- c("p.value")
write.csv(KEGG_enzyme_functions.wilcox1, file = "KEGG_enzyme_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_enzyme_functions2 <- KEGG_enzyme1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_enzyme_functions2 <-  data.frame(ID = 1:ncol(KEGG_enzyme_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_enzyme_functions2$group)
  return(res$p.value)
}

KEGG_enzyme_functions.wilcox2 <-apply(as.matrix(KEGG_enzyme_functions2), 1, fx)
KEGG_enzyme_functions.wilcox2 <- data.frame(KEGG_enzyme_functions.wilcox2 ) 
colnames(KEGG_enzyme_functions.wilcox2) <- c("p.value")
write.csv(KEGG_enzyme_functions.wilcox2, file = "KEGG_enzyme_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_enzyme_functions3 <-KEGG_enzyme1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_enzyme_functions3 <-  data.frame(ID = 1:ncol(KEGG_enzyme_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_enzyme_functions3$group)
  return(res$p.value)
}
KEGG_enzyme_functions.wilcox3 <-apply(as.matrix(KEGG_enzyme_functions3), 1, fx)
KEGG_enzyme_functions.wilcox3 <- data.frame(KEGG_enzyme_functions.wilcox3 ) 
colnames(KEGG_enzyme_functions.wilcox3) <- c("p.value")
write.csv(KEGG_enzyme_functions.wilcox3, file = "KEGG_enzyme_functions.wilcox0.05_T4_T10.csv")

```

```{r}
KEGG_enzyme2=t(KEGG_enzyme1)
KEGG_enzyme3=cbind(KEGG_enzyme2, Colvec_expedition)
pairwise.wilcox.test(KEGG_enzyme3$`ATP-dependent Lon protease [EC:3.4.21.53]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`dipeptidase A [EC:3.4.-.-]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`prolyl oligopeptidase [EC:3.4.21.26]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`protease I [EC:3.2.-.-]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`pyroglutamyl-peptidase [EC:3.4.19.3]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`sortase B`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`dipeptidyl-peptidase 4 [EC:3.4.14.5]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`lactocepin [EC:3.4.21.96]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`two-component system, OmpR family, sensor kinase [EC:2.7.13.3]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`two-component system, AgrA family, sensor histidine kinase AgrC`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`two-component system, NtrC family, sensor histidine kinase HydH`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_enzyme3$`serine/threonine-protein kinase WNK1 [EC:2.7.11.1]`, KEGG_lipid3$Zones,p.adjust.method = "fdr")
```
##enzyme metabolism
```{r , fig.height = 15, fig.width = 10, fig.align = "center"}
KEGG_enzyme_new=read.csv("mean_KEGGFUN_FDR_0.05_enzyme.families.csv",header=T,sep=",")
rownames(KEGG_enzyme_new)=KEGG_enzyme_new[,1]
KEGG_enzyme_new1=KEGG_enzyme_new[,-(1:2)]

KEGG_enzyme_new1_norm <- t(apply(KEGG_enzyme_new1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_enzyme_new1_norm)
pdf(file="Enzyme metabolism",height = 15,width=10)
pheatmap(KEGG_enzyme_new1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE,border_color = NA,cluster_cols = F,cellwidth = 12, cellheight = 15, main="Enzyme Metabolism")

```



```{r}
##cofactors metabolism
KEGG_cofactors=read.csv("mean_KEGGFUN_wilcox_0.05_Metabolism.of.Cofactors.csv",header=T,sep=",")
rownames(KEGG_cofactors)=KEGG_cofactors[,1]
KEGG_cofactors1=KEGG_cofactors[,-(1:2)]
```

```{r}
#Lam and T4
KEGG_cofactors_functions1 <- KEGG_cofactors1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_cofactors_functions1 <-  data.frame(ID = 1:ncol(KEGG_cofactors_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_cofactors_functions1$group)
  return(res$p.value)
}

KEGG_cofactors_functions.wilcox1 <-apply(as.matrix(KEGG_cofactors_functions1), 1, fx)
KEGG_cofactors_functions.wilcox1 <- data.frame(KEGG_cofactors_functions.wilcox1 ) 
colnames(KEGG_cofactors_functions.wilcox1) <- c("p.value")
write.csv(KEGG_cofactors_functions.wilcox1, file = "KEGG_cofactors_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_cofactors_functions2 <- KEGG_cofactors1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_cofactors_functions2 <-  data.frame(ID = 1:ncol(KEGG_cofactors_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_cofactors_functions2$group)
  return(res$p.value)
}

KEGG_cofactors_functions.wilcox2 <-apply(as.matrix(KEGG_cofactors_functions2), 1, fx)
KEGG_cofactors_functions.wilcox2 <- data.frame(KEGG_cofactors_functions.wilcox2 ) 
colnames(KEGG_cofactors_functions.wilcox2) <- c("p.value")
write.csv(KEGG_cofactors_functions.wilcox2, file = "KEGG_cofactors_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_cofactors_functions3 <-KEGG_cofactors1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_cofactors_functions3 <-  data.frame(ID = 1:ncol(KEGG_cofactors_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_cofactors_functions3$group)
  return(res$p.value)
}
KEGG_cofactors_functions.wilcox3 <-apply(as.matrix(KEGG_cofactors_functions3), 1, fx)
KEGG_cofactors_functions.wilcox3 <- data.frame(KEGG_cofactors_functions.wilcox3 ) 
colnames(KEGG_cofactors_functions.wilcox3) <- c("p.value")
write.csv(KEGG_cofactors_functions.wilcox3, file = "KEGG_cofactors_functions.wilcox0.05_T4_T10.csv")
```

```{r}
KEGG_cofactors2=t(KEGG_cofactors1)
KEGG_cofactors3=cbind(KEGG_cofactors2, Colvec_expedition)
pairwise.wilcox.test(KEGG_cofactors3$`dihydroneopterin aldolase [EC:4.1.2.25]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`para-aminobenzoate synthetase component II [EC:2.6.1.85]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`dihydrofolate reductase [EC:1.5.1.3]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`lipoate-protein ligase A [EC:2.7.7.63]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`dephospho-CoA kinase [EC:2.7.1.24]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`type III pantothenate kinase [EC:2.7.1.33]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`adenosylcobinamide-phosphate synthase CobD [EC:6.3.1.10]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`anaerobic magnesium-protoporphyrin IX monomethyl ester cyclase`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`cobalamin biosynthesis protein CbiG`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`acid phosphatase (class A) [EC:3.1.3.2]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`hydroxyethylthiazole kinase [EC:2.7.1.50]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`2-succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylate`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`O-succinylbenzoate synthase [EC:4.2.1.113]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`O-succinylbenzoic acid--CoA ligase [EC:6.2.1.26]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`pyridoxine 4-dehydrogenase [EC:1.1.1.65]`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_cofactors3$`starvation sensing protein RspA`, KEGG_cofactors3$Zones,p.adjust.method = "fdr")

```

##cofactors metabolism
```{r , fig.height = 15, fig.width = 10, fig.align = "center"}
##cofactors metabolism
KEGG_cofactors_new=read.csv("mean_KEGGFUN_FDR_0.05_Metabolism.of.Cofactors.csv",header=T,sep=",")
rownames(KEGG_cofactors_new)=KEGG_cofactors_new[,1]
KEGG_cofactors_new1=KEGG_cofactors_new[,-(1:2)]

KEGG_cofactors_new1_norm <- t(apply(KEGG_cofactors_new1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_cofactors_new1_norm)
pdf(file="Cofactors Metabolism",height = 15,width = 10)
pheatmap(KEGG_cofactors_new1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE,cluster_cols = F,border_color = NA,cellwidth = 12, cellheight = 15,main="Cofactors Metabolism")
```


```{r}
##glycan biosynthesis and metabolism
KEGG_Glycan=read.csv("mean_KEGGFUN_wilcox_0.05_ Glycan.Biosynthesis.csv",header=T,sep=",")
rownames(KEGG_Glycan)=KEGG_Glycan[,1]
KEGG_Glycan1=KEGG_Glycan[,-(1:2)]

KEGG_Glycan1_norm <- t(apply(KEGG_Glycan1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_Glycan1_norm)
pheatmap(KEGG_Glycan1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE,cluster_cols = F,border_color = NA,cellwidth = 12, cellheight = 15,main="Glycan.Biosynthesis and metabolism")
```

```{r}
#Lam and T4
KEGG_Glycan_functions1 <- KEGG_Glycan1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_Glycan_functions1 <-  data.frame(ID = 1:ncol(KEGG_Glycan_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_Glycan_functions1$group)
  return(res$p.value)
}

KEGG_Glycan_functions.wilcox1 <-apply(as.matrix(KEGG_Glycan_functions1), 1, fx)
KEGG_Glycan_functions.wilcox1 <- data.frame(KEGG_Glycan_functions.wilcox1 ) 
colnames(KEGG_Glycan_functions.wilcox1) <- c("p.value")
write.csv(KEGG_Glycan_functions.wilcox1, file = "KEGG_Glycan_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_Glycan_functions2 <- KEGG_Glycan1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_Glycan_functions2 <-  data.frame(ID = 1:ncol(KEGG_Glycan_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_Glycan_functions2$group)
  return(res$p.value)
}

KEGG_Glycan_functions.wilcox2 <-apply(as.matrix(KEGG_Glycan_functions2), 1, fx)
KEGG_Glycan_functions.wilcox2 <- data.frame(KEGG_Glycan_functions.wilcox2 ) 
colnames(KEGG_Glycan_functions.wilcox2) <- c("p.value")
write.csv(KEGG_Glycan_functions.wilcox2, file = "KEGG_Glycan_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_Glycan_functions3 <-KEGG_Glycan1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_Glycan_functions3 <-  data.frame(ID = 1:ncol(KEGG_Glycan_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_Glycan_functions3$group)
  return(res$p.value)
}
KEGG_Glycan_functions.wilcox3 <-apply(as.matrix(KEGG_Glycan_functions3), 1, fx)
KEGG_Glycan_functions.wilcox3 <- data.frame(KEGG_Glycan_functions.wilcox3 ) 
colnames(KEGG_Glycan_functions.wilcox3) <- c("p.value")
write.csv(KEGG_Glycan_functions.wilcox3, file = "KEGG_Glycan_functions.wilcox0.05_T4_T10.csv")
```

```{r}
KEGG_Glycan2=t(KEGG_Glycan1)
KEGG_Glycan3=cbind(KEGG_Glycan2, Colvec_expedition)
pairwise.wilcox.test(KEGG_Glycan3$`poly(glycerol-phosphate) alpha-glucosyltransferase [EC:2.4.1.52]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`UDP-N-acetyl-D-mannosaminuronic acid transferase [EC:2.4.1.-]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`dolichyl-phosphate-mannose-protein mannosyltransferase`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`phosphoheptose isomerase [EC:5.-.-.-]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`galacturonosyltransferase [EC:2.4.1.-]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`undecaprenyl-phosphate galactose phosphotransferase [EC:2.7.8.6]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`alpha-mannosidase [EC:3.2.1.24]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`serine/alanine adding enzyme [EC:2.3.2.10]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_Glycan3$`(heptosyl)LPS beta-1,4-glucosyltransferase [EC:2.4.1.-]`, KEGG_Glycan3$Zones,p.adjust.method = "fdr")

```

##glycan biosynthesis and metabolism

``````{r , fig.height = 15, fig.width = 10, fig.align = "center"}
##glycan biosynthesis and metabolism
KEGG_Glycan_new=read.csv("mean_KEGGFUN_FDR_0.05_Glycan.Biosynthesis.csv",header=T,sep=",")
rownames(KEGG_Glycan_new)=KEGG_Glycan_new[,1]
KEGG_Glycan_new1=KEGG_Glycan_new[,-(1:2)]

KEGG_Glycan_new1_norm <- t(apply(KEGG_Glycan_new1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_Glycan_new1_norm)
pdf(file="Glycan.Biosynthesis and metabolism", height = 15,width=10)
pheatmap(KEGG_Glycan_new1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE,cluster_cols = F,border_color = NA,cellwidth = 12, cellheight = 15,main="Glycan.Biosynthesis and metabolism")
```




```{r}
##other metabolisms
KEGG_others=read.csv("mean_KEGGFUN_wilcox_0.05_Metabolism_others.csv",header=T,sep=",")
rownames(KEGG_others)=KEGG_others[,1]
KEGG_others1=KEGG_others[,-(1:2)]

KEGG_others1_norm <- t(apply(KEGG_others1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_others1_norm)
pheatmap(KEGG_others1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE,border_color = NA,cluster_cols = F,cellwidth = 12, cellheight = 15,main="Other Metabolisms")

```

```{r}
#Lam and T4
KEGG_others_functions1 <- KEGG_others1[,c(idx_Lam, idx_T4)]
clin.sub_KEGG_others_functions1 <-  data.frame(ID = 1:ncol(KEGG_others_functions1), group = rep(c("Lam","T4"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_others_functions1$group)
  return(res$p.value)
}

KEGG_others_functions.wilcox1 <-apply(as.matrix(KEGG_others_functions1), 1, fx)
KEGG_others_functions.wilcox1 <- data.frame(KEGG_others_functions.wilcox1 ) 
colnames(KEGG_others_functions.wilcox1) <- c("p.value")
write.csv(KEGG_others_functions.wilcox1, file = "KEGG_others_functions.wilcox0.05_Lam_T4.csv")


#Lam and T10
KEGG_others_functions2 <- KEGG_others1[,c(idx_Lam, idx_T10)]
clin.sub_KEGG_others_functions2 <-  data.frame(ID = 1:ncol(KEGG_others_functions2), group = rep(c("Lam","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_others_functions2$group)
  return(res$p.value)
}

KEGG_others_functions.wilcox2 <-apply(as.matrix(KEGG_others_functions2), 1, fx)
KEGG_others_functions.wilcox2 <- data.frame(KEGG_others_functions.wilcox2 ) 
colnames(KEGG_others_functions.wilcox2) <- c("p.value")
write.csv(KEGG_others_functions.wilcox2, file = "KEGG_others_functions.wilcox0.05_Lam_T10.csv")

#T10 and T4
KEGG_others_functions3 <-KEGG_Glycan1[,c(idx_T4, idx_T10)]
clin.sub_KEGG_others_functions3 <-  data.frame(ID = 1:ncol(KEGG_others_functions3), group = rep(c("T4","T10"), each=5))

fx <- function(x) {
  res <-wilcox.test(x ~clin.sub_KEGG_others_functions3$group)
  return(res$p.value)
}
KEGG_others_functions.wilcox3 <-apply(as.matrix(KEGG_others_functions3), 1, fx)
KEGG_others_functions.wilcox3 <- data.frame(KEGG_others_functions.wilcox3 ) 
colnames(KEGG_others_functions.wilcox3) <- c("p.value")
write.csv(KEGG_others_functions.wilcox3, file = "KEGG_others_functions.wilcox0.05_T4_T10.csv")

```

```{r}
KEGG_others2=t(KEGG_others1)
KEGG_others3=cbind(KEGG_others2, Colvec_expedition)
pairwise.wilcox.test(KEGG_others3$`2,5-diketo-D-gluconate reductase A [EC:1.1.1.274]`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`carboxylesterase [EC:3.1.1.1]`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`FMN-dependent NADH-azoreductase [EC:1.7.-.-]`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`N-acetylglucosaminyldiphosphoundecaprenol [EC:2.4.1.187]`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`putative acetyltransferase [EC:2.3.1.-]`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`putative colanic acid biosynthesis acetyltransferase WcaF`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`7-alpha-hydroxysteroid dehydrogenase [EC:1.1.1.159]`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`phenylacetic acid degradation protein`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`4-carboxymuconolactone decarboxylase [EC:4.1.1.44]`, KEGG_others3$Zones,p.adjust.method = "fdr")
pairwise.wilcox.test(KEGG_others3$`4-hydroxy-3-methylbut-2-enyl diphosphate reductase [EC:1.17.1.2]`, KEGG_others3$Zones,p.adjust.method = "fdr")
```


```{r , fig.height = 15, fig.width = 10, fig.align = "center"}
##other metabolisms
KEGG_others_new=read.csv("mean_KEGGFUN_FDR_0.05_Metabolism_others.csv",header=T,sep=",")
rownames(KEGG_others_new)=KEGG_others_new[,1]
KEGG_others_new1=KEGG_others_new[,-(1:2)]

KEGG_others_new1_norm <- t(apply(KEGG_others_new1, 1, cal_z_score))
my_sample_col <-  data.frame(sample = rep(c("Lam","T10","T4"), each=5))
row.names(my_sample_col) <- colnames(KEGG_others_new1_norm)
pdf(file="Other Metabolisms",height =15,width=10)
pheatmap(KEGG_others_new1_norm,annotation_col=my_sample_col,show_colnames = F, show_rownames = T,legend = TRUE,border_color = NA,cluster_cols = F,cellwidth = 12, cellheight = 15,main="Other Metabolisms")

```




#Figure S5. Caudovirales phage abundances increased after spinal cord injury. 
```{r}
##Fig 1 a Caudovirales phage abundances are altered after spinal cord injury
Caudovirales_viral_populations=read.csv("Caudovirales_viral_populations.csv",header=T,row.names=NULL,sep=",")
pOTUs2=pOTUs
pOTUs2$viral_contigs=rownames(pOTUs2)

Caudovirales_viral_populations_RPKM=merge(Caudovirales_viral_populations,pOTUs2,by="viral_contigs")

Caudovirales_order=Caudovirales_viral_populations_RPKM[,-(2:3)]
rownames(Caudovirales_order)=Caudovirales_order[,1]
Caudovirales_order=Caudovirales_order[ ,-1]
Caudovirales_order1=data.frame(colSums(Caudovirales_order))
colnames(Caudovirales_order1)=c("Caudovirales")
Caudovirales_order2=cbind(Caudovirales_order1,Colvec_expedition)

#
pairwise.wilcox.test(Caudovirales_order2$Caudovirales, Caudovirales_order2$Zones,p.adjust.method = "fdr")

#Reads mapped to Caudovirales are higher in T4
Caudovirales_order3=ggplot(Caudovirales_order2, aes(x=Zones, y=Caudovirales)) + geom_boxplot(fill=c("#999999", "#E69F00", "#56B4E9"))+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
  labs(title = "Caudovirales",x=NULL,y = "per Million mapped reads (RPKM)")+theme_classic()

Caudovirales_boxplot=Caudovirales_order3+theme(plot.title = element_text(size = 28,hjust=0.5))+theme(axis.text.x = element_text(color="black",size = 20))+theme(axis.text.y = element_text(color="black",size = 18),axis.title.y=element_text(size=20))
Caudovirales_boxplot+scale_y_continuous(expand = c(0,0), limits = c(0,115))+geom_segment(y=105,yend=105,x=1,xend=3)+geom_text(x=2,y=105,label="*",size=10)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

```

#Figure S8. Heat map showing the differential abundance analysis of phages grouped by infected bacterial hosts. 
```{r}
predicted_pOTUs_genus=read.csv("predicted_hosts_genus_RPKM.csv",header=T,sep=",")
rownames(predicted_pOTUs_genus)=predicted_pOTUs_genus[,1]
predicted_pOTUs_genus1=predicted_pOTUs_genus[,-1]

predicted_pOTUs_genus_norm <- t(apply(predicted_pOTUs_genus1, 1, cal_z_score))

pheatmap(predicted_pOTUs_genus_norm,annotation_col = my_sample_col, show_colnames =T,legend = TRUE,cluster_cols=FALSE)
```
