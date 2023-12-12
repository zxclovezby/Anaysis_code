
rm(list = ls())
library(Seurat)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(future)
availableCores()#总核
nbrOfWorkers()#查看使用的核
#调出画板
#### Draw a statistical graph of the number of genes/count number/proportion of mitochondrial genes
library(paletteer)
palettes_d_names
pal=paletteer_d("ggsci::default_igv")[-1]

setwd("")
# Read in data
#read.delim和read.delim2是专门读tab文件的
exp=read.delim("GSE162183_Raw_gene_counts_matrix.tab\\GSE162183_Raw_gene_counts_matrix.tab",header = T,row.names = 1)
exp=as.data.frame(t(exp))
#构造group文件，这个文件至少需要分组（Group）信息，
#在这个分析中我们放入了Sample_ID，Patient和Group。
group= exp[,c(1,2)]
group$Patient <- rownames(group)
group$Patient <- sapply(group$Patient, function(x)unlist(strsplit(x,'_'))[1])
group$Group <- substr(group$Patient,1,4)
group <- group[,-c(1,2)]
group$Sample_ID <- rownames(group)
#exp文件要求行是基因名，列是样本名
exp[1:5,1:5]
exp <- as.data.frame(t(exp))
#构造meta文件，其实就是把上面的group文件换个名，列名也换一下。
sce.meta <- data.frame(Sample_ID=group$Sample_ID,Disease=group$Group,Patient=group$Patient,
                       row.names = colnames(exp))

#构造Seurat对象
scRNA = CreateSeuratObject(counts=exp,
                           meta.data = sce.meta,
                           min.cells = 3, 
                           min.features = 200)
#min.features 一个细胞至少表达200个基因
#min.cells 一个基因至少在3个细胞内表达

seurat_data <- scRNA
#1.计算线粒体基因百分比
#### Mitochondrial gene ratio
seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")

#nFeature_RNA是一个细胞内表达的基因数目
#nCount_RNA是count表达值
VlnPlot(seurat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,cols = pal)
ggsave("1-threePlot.pdf",width = 8,height = 5)
plot1 <- FeatureScatter(seurat_data, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = pal)
ploT1 <- FeatureScatter(seurat_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols = pal)
plot1 + ploT1
ggsave("2-cor.pdf",width = 10,height = 5)

ggdensity(seurat_data@meta.data, x = "nCount_RNA", title = "seurat_data",fill = "#B3D3EC",color = "#B3D3EC",rug = T)
ggsave("3-nCount_RNA.pdf",width = 8,height = 5)
ggdensity(seurat_data@meta.data, x = "nFeature_RNA", title = "seurat_data",fill = "#B3D3EC",color = "#B3D3EC",rug = T)
ggsave("4-nFeature_RNA.pdf",width = 8,height = 5)
ggdensity(seurat_data@meta.data, x = "percent.mt", title = "seurat_data",fill = "#B3D3EC",color = "#B3D3EC",rug = T)
ggsave("5-percent.mt.pdf",width = 8,height = 5)
dev.off()

######################### Detect the resolution parameters of each sample cluster. After the parameters are determined, you can block them without executing [test]
set.resolutions <- seq(0.5, 2, by = 0.1)
pdf(file = "6-PCA-test.pdf")

#### seurat_data filt
seurat_data.pro <- subset(seurat_data, subset = nFeature_RNA > 200 & percent.mt < 30 & nCount_RNA > 500 & nFeature_RNA < 6000) 
plan("multisession", workers = 1) 
seurat_data.pro <- NormalizeData(seurat_data.pro, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_data.pro <- FindVariableFeatures(seurat_data.pro, selection.method = "vst", nfeatures = 2000)
seurat_data.pro <- ScaleData(seurat_data.pro, vars.to.regress = c("nCount_RNA", "percent.mt"))

seurat_data.pro <- RunPCA (seurat_data.pro, features = VariableFeatures(object = seurat_data.pro), ndims.print = 1:2)
ElbowPlot(object = seurat_data.pro, ndims = 100)

seurat_data.pro <- FindNeighbors(seurat_data.pro, dims = 1:40)
seurat_data.pro  <- FindClusters(object = seurat_data.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(seurat_data.pro)
seurat_data.pro  <- RunUMAP(seurat_data.pro , dims = 1:40)
seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = seurat_data.pro, reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
saveRDS(seurat_data.pro,file = "seurat_data.pro.rds")

seurat_data.pro <- readRDS("seurat_data.pro.rds")
DimPlot(object = seurat_data.pro, reduction = 'umap', group.by = "Patient",cols = pal)
ggsave("7-Patient.pdf",width = 7,height = 7)
DimPlot(object = seurat_data.pro, reduction = 'umap', group.by = "Disease",cols = pal)
ggsave("8-Disease.pdf",width = 7,height = 7)

#






