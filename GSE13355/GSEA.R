
#
rm(list = ls())
setwd('')
load('GSE13355.rdata')

#只要保留银屑病患者样本进行分析
table(GSE13355_clin$Tissue)
les_id <- rownames(GSE13355_clin[GSE13355_clin$Tissue=="involved skin from cases",])

expr <- GSE13355_expr[,les_id]
MRGPRX2_expr <- expr["MRGPRX2",]
MRGPRX2_Group <- as.data.frame(t(MRGPRX2_expr))
MRGPRX2_Group$Group <- ifelse(MRGPRX2_Group$MRGPRX2>median(MRGPRX2_Group$MRGPRX2),"High","Low")

cluster <- MRGPRX2_Group[order(MRGPRX2_Group$Group),]#H组在前
table(cluster$Group)#H29,L29
data <- expr[,rownames(cluster)]
identical(rownames(cluster),colnames(data))
range(data)

library(limma)
table(cluster$Group)

group_list <- factor(rep(c('high','low'),c(29,29)))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(data)
contrast.matrix <- makeContrasts('high-low',levels = design)
fit <- lmFit(data,design = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)

alldiff=topTable(fit2,coef = 1,n = Inf)
alldiff <- alldiff[order(alldiff$logFC,decreasing = T),]


library(clusterProfiler)
library(org.Hs.eg.db)
df.id <- bitr(rownames(alldiff), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(df.id)

my <- merge(alldiff,df.id,by.x=0,by.y=1) %>% .[,c(1,2,8)]
colnames(my)[1] <- 'id'
my <- my[order(my$logFC,decreasing = T),]
my$logFC <- 2^my$logFC
gene.expr = my$logFC
names(gene.expr) <- my$ENTREZID
head(gene.expr)

kk <- gseKEGG(gene.expr, organism = "hsa")
go <- gseGO(gene.expr,OrgDb =org.Hs.eg.db)

sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
sortgo<-go[order(go$enrichmentScore, decreasing = T),]

write.table(sortkk,"gsea_kk_output.txt",sep = "\t",quote = F,col.names = T,row.names = F)
write.table(sortgo,"gsea_go_output.txt",sep = "\t",quote = F,col.names = T,row.names = F)

library(ggplot2)
library(enrichplot)
library(paletteer) 
#提供了 R 编程环境中提供的数百种其他调色板的组合集合，详情可以查看此包，总有你满意的方案  
palettes_d_names
pal <- paletteer_d("ggsci::default_igv")[-1] #有几群细胞需要标记就选几种颜色

##选择感兴趣的path
#paths <- c("hsa03030","hsa03430","hsa03040","hsa04110",
#           "hsa00071","hsa00220","hsa04740","hsa00120") 
#gseaplot2(kk, paths,subplots=1:2)
#ggsave("GSEA-KEGG.pdf",width = 6,height = 6)

##选择最显著的前四条和后四条
gseaplot2(kk,row.names(sortkk)[c(1:4,as.numeric(nrow(sortkk)-3):nrow(sortkk))],
          color = pal[1:8],
          subplots=1:2)
ggsave("GSEA-KEGG.pdf",width = 6,height = 6)


##选择感兴趣的gOTerms
#gOTerms<- c("GO:0000727", "GO:1903405","GO:1904851","GO:1904867",
#            "GO:0033539", "GO:0006570","GO:0001867","GO:0007597") 
#gseaplot2(go, gOTerms)
#ggsave("GSEA-GO.pdf",width = 6,height = 6)

##选择最显著的前四条和后四条
gseaplot2(go,row.names(sortgo)[c(1:4,as.numeric(nrow(sortgo)-3):nrow(sortgo))],
          color = pal[1:8],
          subplots=1:2)
ggsave("GSEA-GO.pdf",width = 6,height = 6)

#end


