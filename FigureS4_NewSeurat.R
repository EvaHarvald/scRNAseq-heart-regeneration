rm(list=ls())
install.packages('BiocManager')
library(BiocManager)
BiocManager::install('multtest')
install.packages('Seurat')
library(Seurat)

packageVersion("Seurat")

library(ggplot2)
library(cowplot)
library(usethis)
library(devtools)
library(dplyr)




#######Nextseq
#####E16.5.4


dat <- Read10X(data.dir = "")

colnames(dat) = paste0("E16.4n_", colnames(dat) )
E16.4n <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "10X_E16.4")


#####P5 4n

dat2 <- Read10X(data.dir = "")

colnames(dat2) = paste0("P5.4n_", colnames(dat2) )
P5.4n <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "10X_P5.4")





#################################################
#Combine data from E16.5 and P5:


#Merged Seurat objects.
Merged_seurat_object <- merge(x = E16.4n, y = P5.4n, add.cell.ids = c("E16.4n", "P5.4n"), project = "mouse_combined")
table(Merged_seurat_object$orig.ident)


##-----------------------------------------------------Filtering cells based on markers----------------------------------------------------##
#remove cells not expressing myl7, myh7l, cmlc1, tnnt2a and cells that express >10 myh6.
Merged_seurat_object <- subset(Merged_seurat_object, subset =  Tnni3 > 1 & Tnnt2 > 1 & Tnnc1 > 1 & Actc1 > 1)
table(Merged_seurat_object$orig.ident)


##-----------------------------------------------------Perform Integration----------------------------------------------------##
list <- SplitObject(Merged_seurat_object, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

#Perform the integration
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features, dims = 1:20)
combined <- IntegrateData(anchors, dims = 1:20 )

combined <- ScaleData(combined, vars.to.regress = "nCount_RNA", verbose = TRUE)

DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering

combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.15)


DimPlot(combined, group.by = "orig.ident")


tiff(file = "Umap_clusters_nextseq.tiff", width = 4000, height = 3000, units = "px", res = 600, pointsize = 12)
DimPlot(combined, group.by = "seurat_clusters",  split.by = "orig.ident")
dev.off()

DimPlot(combined, reduction = "umap", split.by = "orig.ident")

combined <- SetIdent(combined, value = "orig.ident" )
Idents(combined)
DefaultAssay(combined) <- "RNA"


#Read in a list of cell cycle markers

cc.genes <- readLines("")

#Aggregate the list into markers for G2/M and S phases
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
DefaultAssay(combined) <- "RNA"

combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, 
                                 set.ident = TRUE)
table(combined$Phase)





Idents(combined) <- "Phase"


combined.G2M <- subset(x = combined, idents =  "G2M")

table(combined.G2M$Phase)



Idents(combined.G2M) <- "orig.ident"

markers <- FindMarkers(combined.G2M, ident.1 = "E16.4n", ident.2 = "P5.4n")
head(markers)
write.csv2(markers, file = "Markers_NextSeq_G2M.csv")



genes <- as.data.frame(combined.G2M@assays[["RNA"]]@meta.features)
write.csv2(genes, file = "Background_NextSeq.csv")




#oPossum results:
#
library(ggplot2)

dt <- read.delim("oPossum_NextSeq_Up in P5.txt", header=T)

dt <- dt[dt$Fisher.score > 0.5,]
dt$TF <- factor(dt$TF, levels = dt$TF[order(dt$Fisher.score)])
ggplot(dt, aes(x= GC.Content, y = Fisher.score)) + 
  geom_point()


dt$TF <- factor(dt$TF, levels = dt$TF[(order(dt$Fisher.score) )])



#Figure S4b:
tiff(file = "oPossum_NextSeq_Up in P5.tiff", width = 3000, height = 4000, units = "px", res = 600, pointsize = 12)
ggplot(dt, aes(y=TF, x=Fisher.score)) + 
  geom_point(stat="identity")+
  geom_point(aes(colour = Target.gene.hits)) +
  theme(legend.position="right")+
  guides(size = guide_legend(order = 1))+
  theme_bw()
dev.off()




dt <- read.delim("oPossum_NextSeq_Up in E16.5.txt", header=T)

dt <- dt[dt$Fisher.score > 0.5,]
dt$TF <- factor(dt$TF, levels = dt$TF[order(dt$Fisher.score)])
ggplot(dt, aes(x= GC.Content, y = Fisher.score)) + 
  geom_point()


dt$TF <- factor(dt$TF, levels = dt$TF[(order(dt$Fisher.score) )])


#Figure S4c:
tiff(file = "oPossum_NextSeq_Up in E16.5.tiff", width = 3000, height = 4000, units = "px", res = 600, pointsize = 12)
ggplot(dt, aes(y=TF, x=Fisher.score)) + 
  geom_point(stat="identity")+
  geom_point(aes(colour = Target.gene.hits)) +
  theme(legend.position="right")+
  guides(size = guide_legend(order = 1))+
  theme_bw()
dev.off()








Cluster_genes_up <- subset(markers,  markers$avg_log2FC>0.2)
Cluster_genes_up <- row.names(Cluster_genes_up)


#### GO Analysis with Clusterprofiler
background <- row.names(combined)

GO_Cluster <- enrichGO(gene  = Cluster_genes_up,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       universe = background,
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.01)

GO_list <- as.data.frame(GO_Cluster)


jpeg("ccGOs_upinE165_G2M_web.jpeg",width=3000,height=2200, units = "px", res = 300)

cnetplot(GO_Cluster, categorySize="pvalue",  showCategory = c("cell division", "positive regulation of cell cycle", "mitotic nuclear division",
                                              "regulation of cell cycle process", "nuclear division", "chromosome segregation",
                                              "organelle fission", "regulation of nuclear division"))
dev.off()



jpeg("ccGOs_upinE16.5_G2M_rainbow.jpeg",width=6000,height=3200, units = "px", res = 300)
cnetplot(GO_Cluster,  circular = TRUE, colorEdge = TRUE,
         showCategory = c("cell division", "positive regulation of cell cycle", "mitotic nuclear division",
                            "regulation of cell cycle process", "nuclear division", "chromosome segregation",
                            "organelle fission", "regulation of nuclear division"))
dev.off()





















##################################
#NovaSeq



#######Novaseq
#####E16.5.4


dat <- Read10X(data.dir = "")

colnames(dat) = paste0("E16.4n_", colnames(dat) )
E16.4n <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "10X_E16.4")

#####P5 4n

dat2 <- Read10X(data.dir = "")


colnames(dat2) = paste0("P5.4n_", colnames(dat2) )
P5.4n <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "10X_P5.4")





#################################################
#Combine data from E16.5 and P5:


#Merged Seurat objects.
Merged_seurat_object <- merge(x = E16.4n, y = P5.4n, add.cell.ids = c("E16.4n", "P5.4n"), project = "mouse_combined")
table(Merged_seurat_object$orig.ident)


##-----------------------------------------------------Filtering cells based on markers----------------------------------------------------##
#remove cells not expressing myl7, myh7l, cmlc1, tnnt2a and cells that express >10 myh6.
Merged_seurat_object <- subset(Merged_seurat_object, subset =  Tnni3 > 1 & Tnnt2 > 1 & Tnnc1 > 1 & Actc1 > 1)
table(Merged_seurat_object$orig.ident)


##-----------------------------------------------------Perform Integration----------------------------------------------------##
list <- SplitObject(Merged_seurat_object, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

#Perform the integration
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features, dims = 1:20)
combined <- IntegrateData(anchors, dims = 1:20 )

combined <- ScaleData(combined, vars.to.regress = "nCount_RNA", verbose = TRUE)

DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering

combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.15)


DimPlot(combined, group.by = "orig.ident")


tiff(file = "Umap_clusters_novaseq.tiff", width = 6000, height = 3000, units = "px", res = 600, pointsize = 12)
DimPlot(combined, group.by = "seurat_clusters",  split.by = "orig.ident")
dev.off()

DimPlot(combined, reduction = "umap", split.by = "orig.ident")

combined <- SetIdent(combined, value = "orig.ident" )
Idents(combined)
DefaultAssay(combined) <- "RNA"


#Read in a list of cell cycle markers

cc.genes <- readLines("")

#Aggregate the list into markers for G2/M and S phases
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
DefaultAssay(combined) <- "RNA"

combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, 
                             set.ident = TRUE)
table(combined$Phase)





Idents(combined) <- "Phase"


combined.G2M <- subset(x = combined, idents =  "G2M")

table(combined.G2M$Phase)



Idents(combined.G2M) <- "orig.ident"

markers <- FindMarkers(combined.G2M, ident.1 = "E16.4n", ident.2 = "P5.4n")
head(markers)
write.csv2(markers, file = "Markers_NovaSeq_G2M.csv")



genes <- as.data.frame(combined.G2M@assays[["RNA"]]@meta.features)
write.csv2(genes, file = "Background_NovaSeq.csv")












#oPossum results:
#
library(ggplot2)

dt <- read.delim("oPossum_NovaSeq_Up in P5.txt", header=T)

dt <- dt[dt$Fisher.score > 0.5,]
dt$TF <- factor(dt$TF, levels = dt$TF[order(dt$Fisher.score)])
ggplot(dt, aes(x= GC.Content, y = Fisher.score)) + 
  geom_point()


dt$TF <- factor(dt$TF, levels = dt$TF[(order(dt$Fisher.score) )])



#Figure S4b:
tiff(file = "oPossum_NovaSeq_Up in P5.tiff", width = 3000, height = 4000, units = "px", res = 600, pointsize = 12)
ggplot(dt, aes(y=TF, x=Fisher.score)) + 
  geom_point(stat="identity")+
  geom_point(aes(colour = Target.gene.hits)) +
  theme(legend.position="right")+
  guides(size = guide_legend(order = 1))+
  theme_bw()
dev.off()




dt <- read.delim("oPossum_NovaSeq_Up in E16.5.txt", header=T)

dt <- dt[dt$Fisher.score > 0.5,]
dt$TF <- factor(dt$TF, levels = dt$TF[order(dt$Fisher.score)])
ggplot(dt, aes(x= GC.Content, y = Fisher.score)) + 
  geom_point()


dt$TF <- factor(dt$TF, levels = dt$TF[(order(dt$Fisher.score) )])


#Figure S4c:
tiff(file = "oPossum_NovaSeq_Up in E16.5.tiff", width=1400,height=1200, units = "px", res = 300)
ggplot(dt[1:20,], aes(y=TF, x=Fisher.score)) + 
  geom_point(stat="identity")+
  geom_point(aes(colour = Target.gene.hits)) +
  theme(legend.position="right")+
  guides(size = guide_legend(order = 1))+
  theme_bw()
dev.off()













Cluster_genes_up <- subset(markers,  markers$avg_log2FC>0.2)
Cluster_genes_up <- row.names(Cluster_genes_up)


#### GO Analysis with Clusterprofiler
background <- row.names(combined)

GO_Cluster <- enrichGO(gene  = Cluster_genes_up,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'SYMBOL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       universe = background,
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.01)

GO_list <- as.data.frame(GO_Cluster)


jpeg("ccGOs_upinE165_G2M_web_Nova.jpeg",width=3000,height=2200, units = "px", res = 300)

cnetplot(GO_Cluster, categorySize="pvalue", showCategory = c("chromosome condensation", 
                                                             "cytoskeleton-dependent cytokinesis", 
                                                             "cytokinesis", "carbohydrate catabolic process",
                                                             "muscle tissue development", "muscle cell apoptotic process"))
dev.off()



jpeg("ccGOs_upinE16.5_G2M_rainbow_Nova.jpeg",width=6000,height=3200, units = "px", res = 300)
cnetplot(GO_Cluster,  circular = TRUE, colorEdge = TRUE, showCategory = c("chromosome condensation", 
                                                                          "cytoskeleton-dependent cytokinesis", 
                                                                          "cytokinesis", "carbohydrate catabolic process",
                                                                          "muscle tissue development", "muscle cell apoptotic process"))
dev.off()














