
#
rm(list = ls())
setwd("")
library(WGCNA)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

load("GSE14905.rdata")
#只要保留银屑病患者样本进行分析
table(GSE14905_clin$Tissue)
les_id <- rownames(GSE14905_clin[GSE14905_clin$Tissue=="lesional",])
expr <- GSE14905_expr[,les_id]

mydata <- expr
dim(mydata)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)

datExpr1<-datExpr0
gsg = goodSamplesGenes(datExpr1, verbose = 3)
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  #想用哪里切，就把“h = 35000”和“cutHeight = 35000”中的“500”换成你的cutoff
  abline(h = 120, col = "red") 

clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10)
keepSamples = (clust==1)
datExpr = datExpr1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

MRGPRX2_expr <- expr["MRGPRX2",]
MRGPRX2_Group <- as.data.frame(t(MRGPRX2_expr))
MRGPRX2_Group$Group <- ifelse(MRGPRX2_Group$MRGPRX2>median(MRGPRX2_Group$MRGPRX2),2,1)
cluster <- MRGPRX2_Group[order(MRGPRX2_Group$Group),]#H组在前
table(cluster$Group)#H16,L17

TNBC_clin <- cluster
new_TNBC_clin <- TNBC_clin[rownames(datExpr),]
new_TNBC_clin$Group <- as.numeric(new_TNBC_clin$Group)
new_TNBC_clin$MRGPRX2 <- as.numeric(new_TNBC_clin$MRGPRX2)

datExpr = as.data.frame(datExpr)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


dev.off()


net = blockwiseModules(datExpr, power = 10,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)

mergedColors = labels2colors(net$colors)

pdf("2module.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]


text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}


samples=new_TNBC_clin

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)

modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.pdf",width = 5.5, height = 13)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.05, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()

#end


