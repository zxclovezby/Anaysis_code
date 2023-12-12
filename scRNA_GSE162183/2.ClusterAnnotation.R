#' @description: annotated the cell type

library(Seurat)
#library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(openxlsx)
library(ggplot2)
library(future)
availableCores()#总核
nbrOfWorkers()#查看使用的核
plan("multisession", workers = 20) 
options(future.globals.maxSize = 100000 * 1024^2) # set 100G RAM
library(paletteer)
palettes_d_names
pal=paletteer_d("ggsci::default_igv")[-1]

setwd("")
source("ratio.plot.R")

seurat_data.pro <- readRDS("seurat_data.pro.rds")

library(tidyr)
library(paletteer)
color_name=palettes_d_names
pal=paletteer_d("khroma::smoothrainbow")

#原始文件注释marker：https://doi.org/10.1038/s41419-021-03724-6
#endothelial (End/Endo), majorly expressing CDH5, PECAM1,CD34; 
#epidermis (EpD), majorly expressing KRT10,KRT14, KRT1; 
#immune (IM), majorly expressing PTPRC(CD45), ITGAM, CD3E; 
#mesenchymal (Mes), majorly expressing PDGFRB, LUM, COL1A1; 
#neural crest-like (Schwann/Melanocyte-like, SchM), majorly expressing SOX10, MLANA, S100B

#End/Endo
FeaturePlot(seurat_data.pro,c("CDH5"),cols = pal)
FeaturePlot(seurat_data.pro,c("PECAM1"),cols = pal) 
FeaturePlot(seurat_data.pro,c("CD34"),cols = pal) 

#EpD
FeaturePlot(seurat_data.pro,c("KRT10"),cols = pal) 
FeaturePlot(seurat_data.pro,c("KRT14"),cols = pal) 
FeaturePlot(seurat_data.pro,c("KRT1"),cols = pal)

#IM
FeaturePlot(seurat_data.pro,c("PTPRC"),cols = pal)
FeaturePlot(seurat_data.pro,c("ITGAM"),cols = pal)
FeaturePlot(seurat_data.pro,c("CD3E"),cols = pal)
#Mast cell
FeaturePlot(seurat_data.pro,c("TPSAB1"),cols = pal)
FeaturePlot(seurat_data.pro,c("TPSB2"),cols = pal)


#Mes
FeaturePlot(seurat_data.pro,c("PDGFRB"),cols = pal)
FeaturePlot(seurat_data.pro,c("LUM"),cols = pal)
FeaturePlot(seurat_data.pro,c("COL1A1"),cols = pal)

#Schwann/Melanocyte-like, SchM
FeaturePlot(seurat_data.pro,c("SOX10"),cols = pal)
FeaturePlot(seurat_data.pro,c("MLANA"),cols = pal)
FeaturePlot(seurat_data.pro,c("S100B"),cols = pal)

#tsrget gene
FeaturePlot(seurat_data.pro,c("TET2"),cols = pal)#都有
FeaturePlot(seurat_data.pro,c("MRGPRX2"),cols = pal)#mast cell一点点


# res 1.5 cluster labels
#seurat_data.pro$RNA_snn_res.1.2=factor(seurat_data.pro$RNA_snn_res.1.2,levels=c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21",	"22",	"23",	"24",	"25",	"26",	"27",	"28",	"29",	"30",	"31",	"32",	"33",	"34",	"35",	"36",	"37",	"38",	"39",	"40",	"41",	"42",	"43",	"44",	"45",	"46",	"47",	"48"))
seurat_data.pro=SetIdent(seurat_data.pro,value="RNA_snn_res.1.2")

new.cluster.ids <- c("Epidermis",
                     "Mesenchymal",
                     "Mesenchymal",
                     "Mesenchymal",
                     "Epidermis",
                     "Epidermis",
                     "Endothelial",
                     "Epidermis",
                     "Immune cell",
                     "Mesenchymal",
                     "Endothelial",
                     "Endothelial",
                     "Immune cell",
                     "Mesenchymal",
                     "Mast cell",
                     "Epidermis",
                     "Endothelial",
                     "Epidermis",
                     "Immune cell",
                     "Mesenchymal",
                     "Immune cell",
                     "Mesenchymal",
                     "Mesenchymal",
                     "Endothelial",
                     "Epidermis",
                     "Mesenchymal",
                     "Epidermis",
                     "Mesenchymal",
                     "Epidermis",
                     "Neural crest-like",
                     "Epidermis",
                     "Mesenchymal",
                     "Epidermis",
                     "Mesenchymal",
                     "Epidermis",
                     "Mesenchymal")
data_merge_labeled <- seurat_data.pro
names(new.cluster.ids) <- levels(data_merge_labeled)
data_merge_labeled <- RenameIdents(data_merge_labeled, new.cluster.ids)
pal=paletteer_d("ggsci::default_igv")[-1]
DimPlot(data_merge_labeled,cols = pal)
saveRDS(data_merge_labeled,"data_merge_labeled.rds")

#---------------------作图-----------------------#
#修改作图函数，可以使得vlnplot添加P值
singlecell_gene_test <- function(SerautObj, 
                                 genes.use, 
                                 group.by=NULL, 
                                 assay = "RNA", 
                                 comp = NULL, 
                                 alpha_start = .05, 
                                 Bonferroni = T,
                                 only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
  else{
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname]
      group1_exp <- group1_exp[which(group1_exp>0)] 
      
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      group2_exp <- group2_exp[which(group2_exp>0)] 
      
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
    
  }
  
  if (Bonferroni == T){
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
    
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
    
  }
  
  else{
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
  }
  
  return(dfOUT)
}

#-------TET2---------#
#分两组，看基因表达的分布
new_seurat <- SetIdent(data_merge_labeled,value = "Disease")
A <- singlecell_gene_test(new_seurat, 
                          genes.use = 'TET2',
                          group.by = 'Disease', 
                          comp = c("Ctrl", "Psor"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(new_seurat, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "Disease",
                         features = 'TET2', 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Ctrl", "Psor"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)

pal=paletteer_d("ggsci::default_igv")[-1]
DimPlot(new_seurat,cols = pal)#"Ctrl", "Psor"分布图
pal=paletteer_d("khroma::smoothrainbow")
FeaturePlot(new_seurat,c("TET2"),cols = pal)#TET2分布图

#提出某类细胞画图
EpD_seurat <- subset(data_merge_labeled,idents = "Epidermis")
EpD_seurat2 <- SetIdent(EpD_seurat,value = "Disease")

A <- singlecell_gene_test(EpD_seurat2, 
                          genes.use = 'TET2',
                          group.by = 'Disease', 
                          comp = c("Ctrl", "Psor"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(EpD_seurat2, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "Disease",
                         features = 'TET2', 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Ctrl", "Psor"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)


#-------MRGPRX2---------#
#分两组，看基因表达的分布
new_seurat <- SetIdent(data_merge_labeled,value = "Disease")
A <- singlecell_gene_test(new_seurat, 
                          genes.use = 'MRGPRX2',
                          group.by = 'Disease', 
                          comp = c("Ctrl", "Psor"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(new_seurat, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "Disease",
                         features = 'MRGPRX2', 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Ctrl", "Psor"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)

pal=paletteer_d("ggsci::default_igv")[-1]
DimPlot(new_seurat,cols = pal)#"Ctrl", "Psor"分布图
pal=paletteer_d("khroma::smoothrainbow")
FeaturePlot(new_seurat,c("MRGPRX2"),cols = pal)#MRGPRX2分布图

#提出某类细胞画图
Mast_seurat <- subset(data_merge_labeled,idents = "Mast cell")
Mast_seurat2 <- SetIdent(Mast_seurat,value = "Disease")

A <- singlecell_gene_test(EpD_seurat2, 
                          genes.use = 'MRGPRX2',
                          group.by = 'Disease', 
                          comp = c("Ctrl", "Psor"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(Mast_seurat2, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "Disease",
                         features = 'MRGPRX2', 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Ctrl", "Psor"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

CombinePlots(plots_violins)


#end
