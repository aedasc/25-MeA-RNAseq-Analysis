#安装包
install.packages('Seurat')
install.packages(c('dplyr','patchwork'))
install.packages('BiocManager')
BiocManager::install('scater')
install.packages('ggplot2')

#加载包
library(Seurat)   
library(dplyr)               
library(patchwork)  
library(cowplot)
library(tidyverse)
library(ggplot2)

#导入数据
pbmc.data <- Read10X(data.dir = "F:/MeA 10X/outs/filtered_feature_bc_matrix")
#min.cells = 3
#min.features = 200
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <-NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
head(pbmc)

#提取表达量变化最显著的10个基因，用于PCA
top10 <- head(VariableFeatures(pbmc), 10)
#可视化高变基因
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#归一化处理
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#PCA降维聚类
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)

#可视化umap
pbmc <- FindNeighbors(pbmc, dims = 1:20)     
pbmc <- FindClusters(pbmc, resolution = 0.3) 
pbmc <- RunUMAP(pbmc, dims = 1:20)
p1 <-DimPlot(pbmc, reduction = "umap")
p1 

#可视化T-SNE
pbmc <- RunTSNE(pbmc, dims = 1:20)
p2 <- DimPlot(pbmc, reduction = "tsne")
p2 

#找marker基因
markers <- FindAllMarkers(pbmc, logfc.threshold = 0.25, min.pct = 0.1, 
                          only.pos = TRUE, pbmc = "wilcox")
pbmc.markers<-markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


#自动区分聚类

library(celldex)
library(SingleR)
library(scater)
library(SummarizedExperiment)

#存储对比用小鼠数据库
Mouse.brain <- MouseRNAseqData()
save(Mouse.brain,file="F:/hip 10X mice+egfp/R/Mouse.brain.RData")
yes
test.seu<-pbmc
test.count=as.data.frame(test.seu[["RNA"]]@counts)


load(file="Mouse.brain.RData")
common_mouse <- intersect(rownames(test.count), rownames(Mouse.brain))
Mouse.brain <- Mouse.brain[common_mouse,]
test.count_forhpca <- test.count[common_mouse,]
test.count_forhpca.se <- SummarizedExperiment(assays=list(counts=test.count_forhpca))
test.count_forhpca.se <- logNormCounts(test.count_forhpca.se)

pred.main.mouse <- SingleR(test = test.count_forhpca.se, ref = Mouse.brain, labels = Mouse.brain$label.main)
result_main_mouse <- as.data.frame(pred.main.mouse$labels)
result_main_mouse$CB <- rownames(pred.main.mouse)
colnames(result_main_mouse) <- c('MOUSE_Main', 'CB')

write.table(result_main_mouse, file = "MOUSE_Main.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE) 
head(result_main_mouse)

test.seu@meta.data$CB=rownames(test.seu@meta.data)
test.seu@meta.data=merge(test.seu@meta.data,result_main_mouse,by="CB")
rownames(test.seu@meta.data)=test.seu@meta.data$CB

p5 <- DimPlot(test.seu, reduction = "tsne", group.by = "MOUSE_Main", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p6 <- DimPlot(test.seu, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
fig_tsne <- plot_grid(p6, p5, labels = c('ident','MOUSE_Main'),rel_widths = c(2,3))
ggsave(filename = "tsne4.pdf", plot = fig_tsne, device = 'pdf', width = 36, height = 12, units = 'cm')

#作图
Oligodendrocyte_progenitor_cells=c(0,1,2)
Oligodendrocyte_cells=c(3)
Endothelial_cells=c(4)
Fibroblasts=c(5)
Macrophages=c(6)
Neutrophils=c(7)
Astrocytes=c(8)
Oligodendrocyte_Astrocytes=c(9)
Pericytes=c(10)

current.cluster.ids <- c(Oligodendrocyte_progenitor_cells,
                         Oligodendrocyte_cells,
                         Endothelial_cells,
                         Fibroblasts,
                         Macrophages,
                         Neutrophils,
                         Astrocytes,
                         Oligodendrocyte_Astrocytes,
                         Pericytes)

new.cluster.ids <- c(rep("Oligodendrocyte progenitor cells",length(Oligodendrocyte_progenitor_cells)),
                     rep("Oligodendrocyte cells",length(Oligodendrocyte_cells)),
                     rep("Endothelial cells",length(Endothelial_cells)),
                     rep("Fibroblasts",length(Fibroblasts)),
                     rep("Macrophages",length(Macrophages)),
                     rep("Neutrophils",length(Neutrophils)),
                     rep("Astrocytes",length(Astrocytes)),
                     rep("Oligodendrocyte/Astrocytes",length(Oligodendrocyte_Astrocytes)),
                     rep("Pericytes",length(Pericytes))
)

test.seu@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(test.seu@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

head(test.seu@meta.data)
table(test.seu@meta.data$celltype)

test.seu@meta.data$celltype<-factor(test.seu@meta.data$celltype,level=c("Oligodendrocyte progenitor cells", 
                                                                        "Oligodendrocyte cells",
                                                                        "Oligodendrocyte/Astrocytes",
                                                                        "Astrocytes",
                                                                        "Endothelial cells",
                                                                        "Pericytes",
                                                                        "Fibroblasts",
                                                                        "Macrophages",
                                                                        "Neutrophils"))

plotCB=as.data.frame(test.seu@meta.data%>%filter(seurat_clusters!="5" & seurat_clusters!="6"& seurat_clusters!="7"))[,"CB"]
DimPlot(test.seu, reduction = "tsne",group.by = "celltype", pt.size=0.3,cells = plotCB,)
ggsave("tense.png",units = "cm",width = 18,height = 10,dpi=1000,bg = "white")  

saveRDS(test.seu,file = "test.seu.rds") #保存test.seu

gene <-read.table("F:/hip 10X mice+egfp/R/gene.txt")

gene <- c("Pdgfra","Pebp1","Oxt","Gpr17","Uchl1","Agrp","Cpt1a","Stat3","Lepr","Cck","Nr4a3")
gene_jyj <- c("Kcnj13","Kcnj10","Kcnj8","Kcnj3","Kcnj9","Kcnj16")
test.seu1<-test.seu

test.seu1@meta.data<-filter(test.seu1@meta.data,seurat_clusters!="5" & seurat_clusters!="6"& seurat_clusters!="7")

DotPlot(test.seu,features = gene_jyj,group.by = "celltype")+ 
  scale_size()
ggsave("dotplot.png",units = "cm",width = 25,height = 10,dpi=1000,bg = "white")  
  

  FeaturePlot(test.seu,features = "Pdgfra",reduction = "tsne",pt.size = 1,cells = plotCB)+
  scale_
