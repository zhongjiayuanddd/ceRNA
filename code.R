


rm(list=ls())

setwd("F:/GSE138693_RAW/topic/results/miRNA/miRNA_analysis")
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

lnc_up_down=read.table("select_up_down_lnc_up_down__analysis_mRNA.txt",header=F)[,1]
lnc_up_downSet=bitr(lnc_up_down, 'SYMBOL', "ENTREZID", "org.Hs.eg.db")[, "ENTREZID"]  
lnc_up_downKEGG <- enrichKEGG(lnc_up_downSet, keyType = "kegg",pvalueCutoff=0.05,qvalueCutoff=0.05,pAdjustMethod = "BH",
                        organism = "hsa", use_internal_data=FALSE)  #KEGG富集分析


enichKK=DOSE::setReadable(lnc_up_downKEGG,OrgDb = org.Hs.eg.db,keyType ="ENTREZID")
lnc_up__down_KEGG_summay<-data.frame(enichKK@result)
write.csv(lnc_up__down_KEGG_summay,"lnc_up_down_KEGG_enichKK.csv")





go <- enrichGO(lnc_up_downSet, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,keyType = 'ENTREZID')

go=DOSE::setReadable(go,OrgDb = org.Hs.eg.db,keyType ="ENTREZID")
Gene_GO_summay<-data.frame(go@result)
write.csv(Gene_GO_summay,"lnc_up_down_GO_enich.csv")





#GO
setwd("F:/GSE138693_RAW/topic/results/miRNA/miRNA_analysis")
library(ggplot2)
library(tidyverse)
data<-read.table("select_GO.txt",header = T,sep = "\t")
data$GO_term<- as.character(data$Description)  #先转换成字符串
data$GO_term<-factor(data$GO_term,levels = c(data$GO_term)) #再强制加入因子
COLS <- c("#FD8D62","#66C3A5", "#8DA1CB")
a<-ggplot(data=data, aes(x=GO_term,y=p, fill=Category)) + geom_bar(stat="identity", width=0.8) + coord_flip() +  xlab("GO term") + ylab("-log10(pvalue)") +scale_fill_manual(values = COLS)+ theme_bw()
ggsave(a, file="go_all.pdf", width=9.03, height=5.74)







rm(list = ls())
setwd("F:/GSE138693_RAW/topic/results/miRNA/miRNA_analysis")
library(ggplot2)
library(ggpubr)
library(Cairo) 

#darkpath<-read.csv("pathway.csv",header = T)
pathway<-read.csv("KEGG_enrich_lnc_up_down.csv",header = T)

CairoPDF("KEGG_enrich_lnc_up_down.pdf",width = 6.8,height =5)
p = ggplot(pathway,aes(GeneRatio,Term))
p=p + geom_point()
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=-1*log10(Pvalue)))
pr = pbubble+scale_color_gradient(low="green",high = "red")
pr = pr+labs(color=expression(-log[10](Pvalue)),size="Count",  
             x="GeneRatio",y="Term")
pr + theme_bw()
dev.off()
##






rm(list=ls())

#无聚类
setwd("F:/GSE138693_RAW/topic/results/miRNA/miRNA_analysis")
library(cluster)
library(pheatmap)
rt<-read.csv("good_survival_heatmap.csv",row.names = 1)
#d<-dist(rt,method = "euclidean")
re<-as.matrix(rt)
re=t(re)


#re[re>1]=1

sample_list=c(rep('Control',1),
              rep('Case',1))

annotation_c <- data.frame(sample_list)
rownames(annotation_c) <- colnames(re)

pdf("plot.pdf",
    width = 8,
    height = 20)
pheatmap(re, #表达数据
         cluster_rows = T,#行聚类
         cluster_cols = F,#列聚类
         annotation_col =annotation_c, #样本分类数据
         annotation_legend=TRUE, # 显示样本分类
         show_rownames = T,# 显示行名
         show_colnames = F,# 显示列名
         scale = 'none', #对行标准化
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100) # 热图基准颜色
)
dev.off()



##
rm(list=ls())

setwd("F:/GSE138693_RAW/topic/results/miRNA/miRNA_analysis")
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

Gene=read.table("pathway_mRNA_good_survival.txt",header=F)[,1]
GeneSet=bitr(Gene, 'SYMBOL', "ENTREZID", "org.Hs.eg.db")[, "ENTREZID"]  
Gene_KEGG <- enrichKEGG(GeneSet, keyType = "kegg",pvalueCutoff=0.05,qvalueCutoff=0.05,pAdjustMethod = "BH",
                             organism = "hsa", use_internal_data=FALSE)  #KEGG富集分析

enichKK=DOSE::setReadable(Gene_KEGG,OrgDb = org.Hs.eg.db,keyType ="ENTREZID")

summay<-data.frame(enichKK@result)

flodlist<-read.csv("cnetplot_2flod.csv",header = T)

genename<-flodlist[,1]


flodlist=flodlist[,-1]
flodlist[flodlist<(-20)]=-10


flodlist=as.matrix(flodlist)
rownames(flodlist)<-genename

foldChange<-flodlist
foldChange[,1]

library(Cairo)
library(ggplot2)
CairoPDF("cnetplot.pdf",width = 12,height = 10)
p1<-cnetplot(enichKK,foldChange =foldChange[,1] )
p1+scale_color_gradientn(colours = c("green","#FF8C00", "red"))
dev.off()







rm(list = ls())
setwd("F:/GSE138693_RAW/topic/results/miRNA/miRNA_analysis/ceRNA")
library(Cairo) 
b<-matrix(c(886,1533,580,266),nrow = 2,ncol = 2)
CairoPDF("barplot1.pdf",width =5,height =6)
barplot(b,col=c('#F5F5DC','#FFAEB9'),beside=F,ylim =c(0,3000))
dev.off()


























