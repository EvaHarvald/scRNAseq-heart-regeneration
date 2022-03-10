library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat)
install.packages("magrittr") # only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # need to run every time you start R and want to use %>%
library(dplyr)    # alternative, this also loads %>%
rm(list=ls())


##################### Cluster from E16.5 vs. P5 #########################


#E16.5 2n
dat <- Read10X(data.dir = "")
colnames(dat) = paste0("E16.5.2n_", colnames(dat) )
E165.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")

# Filter out cells with no expression of Tnni3:


E16.5.2n.filter <- subset(E165.2n, subset = Tnni3> 1)
E16.5.2n.filter <- subset(E165.2n, subset = Tnnt2>1)
E16.5.2n.filter <- subset(E165.2n, subset = Tnnc1>1)
E16.5.2n.filter <- subset(E165.2n, subset = Actc1> 1)
table(E16.5.2n.filter@meta.data$orig.ident)


# normalize and identify variable features for each dataset independently

E16.5.2n.filter <- NormalizeData(E16.5.2n.filter)
E16.5.2n.filter <- FindVariableFeatures(E16.5.2n.filter, selection.method = "vst", nfeatures = 2000)

E16.5.2n.filter <- ScaleData(E16.5.2n.filter, vars.to.regress = "nCount_RNA", verbose = TRUE)


# Run the standard workflow for visualization and clustering
E16.5.2n.filter <- RunPCA(E16.5.2n.filter, npcs = 30, verbose = FALSE)
E16.5.2n.filter <- RunUMAP(E16.5.2n.filter, reduction = "pca", dims = 1:20)
E16.5.2n.filter <- FindNeighbors(E16.5.2n.filter, reduction = "pca", dims = 1:20)
E16.5.2n.filter <- FindClusters(E16.5.2n.filter, resolution = 0.15)







#E16.5 4n
dat1 <- Read10X(data.dir = "")
colnames(dat1) = paste0("E16.5.4n_", colnames(dat1) )
E165.4n <- CreateSeuratObject(counts = dat1, min.cells = 5, min.features = 200, project = "E165_4n")

# Filter out cells with no expression of Tnni3:


E16.5.4n.filter <- subset(E165.4n, subset = Tnni3> 1)
E16.5.4n.filter <- subset(E165.4n, subset = Tnnt2>1)
E16.5.4n.filter <- subset(E165.4n, subset = Tnnc1>1)
E16.5.4n.filter <- subset(E165.4n, subset = Actc1> 1)
table(E16.5.4n.filter@meta.data$orig.ident)


# normalize and identify variable features for each dataset independently

E16.5.4n.filter <- NormalizeData(E16.5.4n.filter)
E16.5.4n.filter <- FindVariableFeatures(E16.5.4n.filter, selection.method = "vst", nfeatures = 2000)

E16.5.4n.filter <- ScaleData(E16.5.4n.filter, vars.to.regress = "nCount_RNA", verbose = TRUE)


# Run the standard workflow for visualization and clustering
E16.5.4n.filter <- RunPCA(E16.5.4n.filter, npcs = 30, verbose = FALSE)
E16.5.4n.filter <- RunUMAP(E16.5.4n.filter, reduction = "pca", dims = 1:20)
E16.5.4n.filter <- FindNeighbors(E16.5.4n.filter, reduction = "pca", dims = 1:20)
E16.5.4n.filter <- FindClusters(E16.5.4n.filter, resolution = 0.15)



#Merged Seurat objects.
Merged_seurat_object.E16 <- merge(x = E16.5.2n.filter, y = E16.5.4n.filter, add.cell.ids = c("E16.2n", "E16.4n"), project = "mouse_combined")
table(Merged_seurat_object.E16$orig.ident)



##-----------------------------------------------------Perform Integration----------------------------------------------------##
list <- SplitObject(Merged_seurat_object.E16, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

#Perform the integration
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features, dims = 1:20)
combined.E16 <- IntegrateData(anchorset = anchors, dims = 1:20 )

combined.E16 <- ScaleData(combined.E16, vars.to.regress = "nCount_RNA", verbose = TRUE)

DefaultAssay(combined.E16) <- "integrated"

# Run the standard workflow for visualization and clustering

combined.E16 <- RunPCA(combined.E16, npcs = 30, verbose = FALSE)
combined.E16 <- RunUMAP(combined.E16, reduction = "pca", dims = 1:20)
combined.E16 <- FindNeighbors(combined.E16, reduction = "pca", dims = 1:20)
combined.E16 <- FindClusters(combined.E16, resolution = 0.15)

#Figure b:

UMAPPlot(combined.E16)

#-----------------------------------------------------------------------------------------------------------------------------------
#Read in a list of cell cycle markers

cc.genes <- readLines("")

#Aggregate the list into markers for G2/M and S phases
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
DefaultAssay(combined.E16) <- "RNA"

combined.E16 <- CellCycleScoring(combined.E16, s.features = s.genes, g2m.features = g2m.genes, 
                                set.ident = TRUE)
table(combined.E16$Phase)

##Alter threshold value
#combined.E16$Phase[combined.E16$S.Score < 0.1 & combined.E16$G2M.Score < 0.1] <- "G1"

table(combined.E16$Phase)


Idents(combined.E16) <- "Phase"


p1 <- UMAPPlot(combined.E16,  pt.size = 1)


p2 <- DimPlot(combined.E16, reduction = "umap", split.by = "orig.ident", group.by = "Phase")



prop.table(table(Idents(combined.E16)))

table(Idents(combined.E16), combined.E16$orig.ident)

table(combined.E16$orig.ident)

d <- prop.table(table(Idents(combined.E16), combined.E16$orig.ident), margin = 2)

p3 <- tableGrob(format(round(d, 3)))

grid.arrange(p3)




tiff(file = "_CC_E165_correctsettings.tiff", width = 12000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1, p2, p3, ncol=3 ,widths=c(6,8,2)))

dev.off()
##############################################################




Idents(combined.E16) <- "seurat_clusters"


####################ClusterProfiler###########################

combined.E165_markers <- FindAllMarkers(combined.E16, min.pct = 0.25, log2fc.threshold = 0.25)
head(combined.E165_markers)

write.csv2(combined.E165_markers, file="all_markers_E165.csv")


Num_Clusters <- unique(combined.E165_markers$cluster)

Cluster_gene_list <- list()
for (i in Num_Clusters){ 
  Cluster_genes <- subset(combined.E165_markers, combined.E165_markers$cluster == i & combined.E165_markers$avg_log2FC>0.5)
  Cluster_genes <- Cluster_genes %>% pull(gene)
  Cluster_gene_list[[i]] <- Cluster_genes
}




#### GO Analysis with Clusterprofiler
background <- row.names(combined.E16)

GO_list <- list()
GO_list_simplified <- list()
for (i in Num_Clusters){
  GO_Cluster <- enrichGO(gene  = Cluster_gene_list[[i]],
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         universe = background,
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.01)
  GO_list[[i]] <- as.data.frame(GO_Cluster)
  GO_list_simplified[[i]] <- as.data.frame(simplify(GO_Cluster, cutoff=0.6, by="p.adjust", select_fun=min))
  
}


top_GO <- list()

for (i in Num_Clusters){
  
  top_GO <- lapply(GO_list[["i"]], "[[",  "Description[[1]]")
  
}


cat(capture.output(print(GO_list), file="GO_list_E165.txt"))


top_GO <- list()

top_GO <- lapply(GO_list, `[[`, 2)

top_GO_first <- list()
top_GO_first <- lapply(top_GO, `[`, 4)


new.cluster.ids <- sapply(1:length(top_GO_first), function(i) top_GO_first[[i]])



new.cluster.ids <- make.names(new.cluster.ids,unique=T)


names(new.cluster.ids) <- levels(combined.E16)

combined.E16 <- RenameIdents(combined.E16, new.cluster.ids)


p1 <- UMAPPlot(combined.E16,  pt.size = 1)


p2 <- DimPlot(combined.E16, reduction = "umap", split.by = "orig.ident", group.by = "seurat_clusters")





prop.table(table(Idents(combined.E16)))

table(Idents(combined.E16), combined.E16$orig.ident)

table(combined.E16$orig.ident)

d <- prop.table(table(Idents(combined.E16), combined.E16$orig.ident), margin = 2)

p3 <- tableGrob(format(round(d, 2)))

grid.arrange(p3)




tiff(file = "_GOterms_E165.tiff", width = 13000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1 + NoLegend(), p2, p3, ncol=3 ,widths=c(6,10,6)))

dev.off()


##########################



tiff(file = "_sample_E165.tiff", width = 12000, height = 4000, units = "px", res = 650, pointsize = 12)

DimPlot(combined.E16, reduction = "umap", label = TRUE, split.by="orig.ident", label.size = 8,  pt.size = 2)

dev.off()








Idents(combined.E16) <- "seurat_clusters"



t=combined.E16@meta.data
summary(as.factor(t$old.ident)) 

combined.E16.5.subset <- subset(x = combined.E16, idents = c("1","3"))

t=combined.E16.5.subset@meta.data
summary(as.factor(t$orig.ident)) 





combined.E16.5.subset@meta.data$Phase <- paste0(combined.E16.5.subset@meta.data$Phase)

Idents(combined.E16.5.subset) <- "Phase"

t=combined.E16.5.subset@meta.data
summary(t$Phase)

combined.subset.G2M.E16.5 <- subset(x = combined.E16.5.subset, idents =  "G2M")
combined.subset.S.E16.5 <- subset(x = combined.E16.5.subset, idents =  "S")


table(combined.subset.G2M.E16.5@meta.data$orig.ident)
#E16.5.2n E16.5.4n 
#387      1778 
table(combined.subset.S.E16.5@meta.data$orig.ident)
#E16.5.2n E16.5.4n 
#959     297 



cells.E165 = rownames(combined.subset.G2M.E16.5@meta.data[combined.subset.G2M.E16.5@meta.data$orig.ident== "E16.5.4n",])


tiff(file = "_highlighted_cells.E165.tiff", width = 5000, height = 4000, units = "px", res = 650, pointsize = 12)

UMAPPlot(combined.E16,cells.highlight=cells.E165)
dev.off()






Idents(combined.E16) <- "seurat_clusters"


cells.E165 = rownames(combined.E16@meta.data[combined.E16@meta.data$seurat_clusters==c("1","3"),])


tiff(file = "_highlighted_E165.tiff", width = 5000, height = 4000, units = "px", res = 650, pointsize = 12)

UMAPPlot(combined.E16,cells.highlight=cells.E165)
dev.off()

###########################################################################




#P5 2n
dat2 <- Read10X(data.dir = "")
colnames(dat2) = paste0("P5.2n_", colnames(dat2) )
P5.2n  <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "mouse-P5_2n")

# Filter out cells with no expression of Tnni3:


P5.2n.filter <- subset(P5.2n, subset = Tnni3> 1)
P5.2n.filter <- subset(P5.2n, subset = Tnnt2>1)
P5.2n.filter <- subset(P5.2n, subset = Tnnc1>1)
P5.2n.filter <- subset(P5.2n, subset = Actc1> 1)
table(P5.2n.filter@meta.data$orig.ident)


# normalize and identify variable features for each dataset independently

P5.2n.filter <- NormalizeData(P5.2n.filter)
P5.2n.filter <- FindVariableFeatures(P5.2n.filter, selection.method = "vst", nfeatures = 2000)

P5.2n.filter <- ScaleData(P5.2n.filter, vars.to.regress = "nCount_RNA", verbose = TRUE)


# Run the standard workflow for visualization and clustering
P5.2n.filter <- RunPCA(P5.2n.filter, npcs = 30, verbose = FALSE)
P5.2n.filter <- RunUMAP(P5.2n.filter, reduction = "pca", dims = 1:20)
P5.2n.filter <- FindNeighbors(P5.2n.filter, reduction = "pca", dims = 1:20)
P5.2n.filter <- FindClusters(P5.2n.filter, resolution = 0.15)







#P5 4n
dat3 <- Read10X(data.dir = "")
colnames(dat3) = paste0("P5.4n_", colnames(dat3) )
P5.4n <- CreateSeuratObject(counts = dat3, min.cells = 5, min.features = 200, project = "P5_4n")

# Filter out cells with no expression of Tnni3:


P5.4n.filter <- subset(P5.4n, subset = Tnni3> 1)
P5.4n.filter <- subset(P5.4n, subset = Tnnt2>1)
P5.4n.filter <- subset(P5.4n, subset = Tnnc1>1)
P5.4n.filter <- subset(P5.4n, subset = Actc1> 1)
table(P5.4n.filter@meta.data$orig.ident)


# normalize and identify variable features for each dataset independently

P5.4n.filter <- NormalizeData(P5.4n.filter)
P5.4n.filter <- FindVariableFeatures(P5.4n.filter, selection.method = "vst", nfeatures = 2000)

P5.4n.filter <- ScaleData(P5.4n.filter, vars.to.regress = "nCount_RNA", verbose = TRUE)


# Run the standard workflow for visualization and clustering
P5.4n.filter <- RunPCA(P5.4n.filter, npcs = 30, verbose = FALSE)
P5.4n.filter <- RunUMAP(P5.4n.filter, reduction = "pca", dims = 1:20)
P5.4n.filter <- FindNeighbors(P5.4n.filter, reduction = "pca", dims = 1:20)
P5.4n.filter <- FindClusters(P5.4n.filter, resolution = 0.15)



#Merged Seurat objects.
Merged_seurat_object.P5 <- merge(x = P5.2n.filter, y = P5.4n.filter, add.cell.ids = c("P5.2n", "P5.4n"), project = "mouse_combined")
table(Merged_seurat_object.P5$orig.ident)



##-----------------------------------------------------Perform Integration----------------------------------------------------##
list <- SplitObject(Merged_seurat_object.P5, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

#Perform the integration
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features, dims = 1:20)
combined.P5 <- IntegrateData(anchorset = anchors, dims = 1:20 )

combined.P5 <- ScaleData(combined.P5, vars.to.regress = "nCount_RNA", verbose = TRUE)

DefaultAssay(combined.P5) <- "integrated"

# Run the standard workflow for visualization and clustering

combined.P5 <- RunPCA(combined.P5, npcs = 30, verbose = FALSE)
combined.P5 <- RunUMAP(combined.P5, reduction = "pca", dims = 1:20)
combined.P5 <- FindNeighbors(combined.P5, reduction = "pca", dims = 1:20)
combined.P5 <- FindClusters(combined.P5, resolution = 0.15)




#-----------------------------------------------------------------------------------------------------------------------------------
#Read in a list of cell cycle markers

cc.genes <- readLines("")

#Aggregate the list into markers for G2/M and S phases
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
DefaultAssay(combined.P5) <- "RNA"

combined.P5 <- CellCycleScoring(combined.P5, s.features = s.genes, g2m.features = g2m.genes, 
                             set.ident = TRUE)
table(combined.P5$Phase)

##Alter threshold value
#combined.P5$Phase[combined.P5$S.Score < 0.1 & combined.P5$G2M.Score < 0.1] <- "G1"

table(combined.P5$Phase)


Idents(combined.P5) <- "Phase"


combined.P5@meta.data[["Phase"]] <- factor(combined.P5@meta.data[["Phase"]], levels = c("G1", "S", "G2M"))


p1 <- UMAPPlot(combined.P5,  pt.size = 1)


p2 <- DimPlot(combined.P5, reduction = "umap", split.by = "orig.ident", group.by = "Phase")



prop.table(table(Idents(combined.P5)))

table(Idents(combined.P5), combined.P5$orig.ident)

table(combined.P5$orig.ident)

d <- prop.table(table(Idents(combined.P5), combined.P5$orig.ident), margin = 2)

p3 <- tableGrob(format(round(d, 3)))

grid.arrange(p3)




tiff(file = "_CC_P5_CorrectSettings.tiff", width = 12000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1, p2, p3, ncol=3 ,widths=c(6,8,2)))

dev.off()
##############################################################



Idents(combined.P5) <- "seurat_clusters"

####################ClusterProfiler###########################

combined_markers.P5 <- FindAllMarkers(combined.P5, min.pct = 0.25, log2fc.threshold = 0.25)
head(combined_markers.P5)

write.csv2(combined_markers.P5, file="all_markers_P5.csv")


Num_Clusters <- unique(combined_markers.P5$cluster)

Cluster_gene_list <- list()
for (i in Num_Clusters){ 
  Cluster_genes <- subset(combined_markers.P5, combined_markers.P5$cluster == i & combined_markers.P5$avg_log2FC>0.5)
  Cluster_genes <- Cluster_genes %>% pull(gene)
  Cluster_gene_list[[i]] <- Cluster_genes
}




#### GO Analysis with Clusterprofiler
background <- row.names(combined.P5)

GO_list <- list()
GO_list_simplified <- list()
for (i in Num_Clusters){
  GO_Cluster <- enrichGO(gene  = Cluster_gene_list[[i]],
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         universe = background,
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.01)
  GO_list[[i]] <- as.data.frame(GO_Cluster)
  GO_list_simplified[[i]] <- as.data.frame(simplify(GO_Cluster, cutoff=0.6, by="p.adjust", select_fun=min))
  
}


top_GO <- list()

for (i in Num_Clusters){
  
  top_GO <- lapply(GO_list[["i"]], "[[",  "Description[[1]]")
  
}


cat(capture.output(print(GO_list), file="GO_list_P5.txt"))


top_GO <- list()

top_GO <- lapply(GO_list, `[[`, 2)

top_GO_first <- list()
top_GO_first <- lapply(top_GO, `[`, 1)


new.cluster.ids <- sapply(1:length(top_GO_first), function(i) top_GO_first[[i]])



new.cluster.ids <- make.names(new.cluster.ids,unique=T)


names(new.cluster.ids) <- levels(combined.P5)

combined.P5 <- RenameIdents(combined.P5, new.cluster.ids)


p1 <- UMAPPlot(combined.P5,  pt.size = 1)


p2 <- DimPlot(combined.P5, reduction = "umap", split.by = "orig.ident", group.by = "seurat_clusters")





prop.table(table(Idents(combined.P5)))

table(Idents(combined.P5), combined.P5$orig.ident)

table(combined.P5$orig.ident)

d <- prop.table(table(Idents(combined.P5), combined.P5$orig.ident), margin = 2)

p3 <- tableGrob(format(round(d, 2)))

grid.arrange(p3)



tiff(file = "_GOterms_P5.tiff", width = 13000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1 + NoLegend(), p2, p3, ncol=3 ,widths=c(6,10,6)))

dev.off()

##########################



tiff(file = "_sample_P5.tiff", width = 12000, height = 4000, units = "px", res = 650, pointsize = 12)

DimPlot(combined.P5, reduction = "umap", label = TRUE, split.by="orig.ident", label.size = 8,  pt.size = 2)

dev.off()







Idents(combined.P5) <- "seurat_clusters"




t=combined.P5@meta.data
summary(as.factor(t$orig.ident)) 

combined.P5.subset <- subset(x = combined.P5, idents = c("1", "3"))

t=combined.P5.subset@meta.data
summary(as.factor(t$orig.ident)) 





combined.P5.subset@meta.data$Phase <- paste0(combined.P5.subset@meta.data$Phase)

Idents(combined.P5.subset) <- "Phase"

t=combined.P5.subset@meta.data
summary(t$Phase)

combined.subset.G2M.P5 <- subset(x = combined.P5.subset, idents =  "G2M")
combined.subset.S.P5 <- subset(x = combined.P5.subset, idents =  "S")


table(combined.subset.G2M.P5@meta.data$orig.ident)
#P5.2n P5.4n 
#353      149 
table(combined.subset.S.P5@meta.data$orig.ident)
#P5.2n P5.4n 
#429     66 



cells.P5 = rownames(combined.subset.G2M.P5@meta.data[combined.subset.G2M.P5@meta.data$orig.ident== "P5.4n",])


tiff(file = "_highlighted_P5.tiff", width = 5000, height = 4000, units = "px", res = 650, pointsize = 12)

UMAPPlot(combined.P5,cells.highlight=cells.P5)
dev.off()





########### Merge G2M 

#Merged Seurat objects.
Merged_seurat_object.subset <- merge(x = combined.subset.G2M.E16.5, y = combined.subset.G2M.P5, add.cell.ids = c("E16.5", "P5"), project = "mouse_combined")
table(Merged_seurat_object.subset$orig.ident)



##########################



table(Merged_seurat_object.subset$orig.ident)



Merged_seurat_object.subset <- NormalizeData(Merged_seurat_object.subset)
Merged_seurat_object.subset <- FindVariableFeatures(Merged_seurat_object.subset, selection.method = "vst", nfeatures = 2000)

Merged_seurat_object.subset <- ScaleData(Merged_seurat_object.subset , vars.to.regress = "nCount_RNA", verbose = TRUE)

Merged_seurat_object.subset  <- RunPCA(Merged_seurat_object.subset , npcs = 30, verbose = FALSE)
Merged_seurat_object.subset  <- RunUMAP(Merged_seurat_object.subset , reduction = "pca", dims = 1:20)
Merged_seurat_object.subset  <- FindNeighbors(Merged_seurat_object.subset , reduction = "pca", dims = 1:20)
Merged_seurat_object.subset  <- FindClusters(Merged_seurat_object.subset , resolution = 0.15)

Idents(Merged_seurat_object.subset) <- "orig.ident"


tiff(file = "_selectedInNewUmap.tiff", width = 4000, height = 4000, units = "px", res = 650, pointsize = 12)
UMAPPlot(Merged_seurat_object.subset, pt.size = 2)
dev.off()

table(Merged_seurat_object.subset$orig.ident)

Idents(Merged_seurat_object.subset) <- "Phase"

table(Merged_seurat_object.subset$orig.ident)


######################################

Merged_seurat_object.subset  <- RunUMAP(Merged_seurat_object.subset , reduction = "pca", n.components = 3,  dims = 1:20)

tsne_1 <- Merged_seurat_object.subset@reductions[["umap"]]@cell.embeddings[,1]

tsne_2 <- Merged_seurat_object.subset@reductions[["umap"]]@cell.embeddings[,2]

tsne_3 <- Merged_seurat_object.subset@reductions[["umap"]]@cell.embeddings[,3]


library(scatterplot3d)


scatterplot3d(x = tsne_1, y = tsne_2, z = tsne_3, color = as.numeric(1:4)[Merged_seurat_object.subset@active.ident], pch = 16)


library(rgl) #interactive 3d plotting
library(rmarkdown)
#Figure 4c:


plot3d(x = tsne_1, y = tsne_2, z = tsne_3, col = as.numeric(1:4)[Merged_seurat_object.subset@active.ident], size=1, type='s')

legend3d("topright", legend = levels(Merged_seurat_object.subset@active.ident), col = c("black","red", "green",  "blue"),  pch = 16, cex=1, inset=c(0.02))



##############################################################################################################################


######### For 4n analysis ########


combined.G2M.4n <- subset(x = Merged_seurat_object.subset, idents =  c("E16.5.4n", "P5.4n"))


G2M.response <- FindMarkers(combined.G2M.4n, ident.1 = "E16.5.4n", ident.2 = "P5.4n", 
                            print.bar = FALSE)

head(G2M.response, 15)

write.csv2(G2M.response, file = "Markers between E16.5 and P5_cluster_G2M_wo2n.csv")

table(combined.G2M.4n$orig.ident)

Idents(combined.G2M.4n) <- "orig.ident"

tiff(file = "_selectedInNewUmap_G2M.tiff", width = 4000, height = 4000, units = "px", res = 650, pointsize = 12)
UMAPPlot(combined.G2M.4n, pt.size = 3)
dev.off()




#High in E16.5

gene.use.high.in.E165 = subset(G2M.response,avg_log2FC>0.5&p_val_adj<0.001)

gene.use.high.in.E165.names <- row.names(gene.use.high.in.E165)

DoHeatmap(object = combined.G2M.4n, features = gene.use.high.in.E165.names)




#High in P5


gene.use.high.in.P5 <- subset(G2M.response,avg_log2FC< -0.5 &p_val_adj<0.001)

gene.use.high.in.P5.names <- row.names(gene.use.high.in.P5)

DoHeatmap(object = combined.G2M.4n, features = gene.use.high.in.P5)






#############################################################################################################################


###################################################################################################


############################## GO terms ##########################




#### GO Analysis with Clusterprofiler
background <- row.names(combined.G2M.4n)

GO.terms.up.in.E165 <- enrichGO(gene  = gene.use.high.in.E165.names,
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         universe = background,
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.01)



  
head(GO.terms.up.in.E165@result[["Description"]])
  




jpeg("GO.terms.up.in.E165_web.jpeg",width=3000,height=2200, units = "px", res = 300)

cnetplot(GO.terms.up.in.E165, categorySize="pvalue", foldChange=gene.use.high.in.E165.names)
dev.off()



jpeg("GO.terms.up.in.E165_rainbow.jpeg",width=6000,height=3200, units = "px", res = 300)

cnetplot(GO.terms.up.in.E165, foldChange=gene.use.high.in.E165.names, circular = TRUE, colorEdge = TRUE)
dev.off()



#Up in P5:

GO.terms.up.in.P5 <- enrichGO(gene  = gene.use.high.in.P5.names,
                                OrgDb         = org.Mm.eg.db,
                                keyType       = 'SYMBOL',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                universe = background,
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.01)






head(GO.terms.up.in.P5@result[["Description"]])



jpeg("GO.terms.up.in.P5_web.jpeg",width=3000,height=2200, units = "px", res = 300)

cnetplot(GO.terms.up.in.P5, categorySize="pvalue", foldChange=gene.use.high.in.P5.names)
dev.off()



jpeg("GO.terms.up.in.P5_rainbow.jpeg",width=6000,height=3200, units = "px", res = 300)

cnetplot(GO.terms.up.in.P5, foldChange=gene.use.high.in.P5.names, circular = TRUE, colorEdge = TRUE)
dev.off()





############### TFs POSSUM #################

write.csv2(background, file = "Background.csv")









dt <- read.delim("TFs_upinclusterE165vsP5_possum.txt", header=T)

dt <- dt[dt$Fisher.score > 0.5,]
dt$TF <- factor(dt$TF, levels = dt$TF[order(dt$Fisher.score)])
ggplot(dt, aes(x= GC.Content, y = Fisher.score)) + 
  geom_point()


dt$TF <- factor(dt$TF, levels = dt$TF[(order(dt$Fisher.score) )])

#Figure 4e:
tiff("TFs_upinclusterE165vsP5_possum.tiff",width=1400,height=1200, units = "px", res = 300)
ggplot(dt[1:20,], aes(y=TF, x=Fisher.score)) + 
  geom_point(stat="identity")+
  geom_point(aes(colour = Target.gene.hits)) +
  theme(legend.position="right")+
  guides(size = guide_legend(order = 1))+
  theme_bw()
dev.off()



dt <- read.delim("TFs_upinclusterP5vsE156_possum.txt", header=T)

dt <- dt[dt$Fisher.score > 0.5,]
dt$TF <- factor(dt$TF, levels = dt$TF[order(dt$Fisher.score)])
ggplot(dt, aes(x= GC.Content, y = Fisher.score)) + 
  geom_point()


dt$TF <- factor(dt$TF, levels = dt$TF[(order(dt$Fisher.score) )])


tiff("TFs_upinclusterP5vsE156_possum.tiff",width=2000,height=2200, units = "px", res = 300)
ggplot(dt, aes(y=TF, x=Fisher.score)) + 
  geom_point(stat="identity")+
  geom_point(aes(colour = Target.gene.hits)) +
  theme(legend.position="right")+
  guides(size = guide_legend(order = 1))+
  theme_bw()
dev.off()





############################################ P1 ###############################


#P1 2n
dat2 <- Read10X(data.dir = "")
colnames(dat2) = paste0("P1.2n_", colnames(dat2) )
P1.2n  <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "mouse-P1_2n")

# Filter out cells with no expression of Tnni3:


P1.2n.filter <- subset(P1.2n, subset = Tnni3> 1)
P1.2n.filter <- subset(P1.2n, subset = Tnnt2>1)
P1.2n.filter <- subset(P1.2n, subset = Tnnc1>1)
P1.2n.filter <- subset(P1.2n, subset = Actc1> 1)
table(P1.2n.filter@meta.data$orig.ident)




#P1 4n
dat3 <- Read10X(data.dir = "")
colnames(dat3) = paste0("P1.4n_", colnames(dat3) )
P1.4n <- CreateSeuratObject(counts = dat3, min.cells = 5, min.features = 200, project = "P1_4n")

# Filter out cells with no expression of Tnni3:


P1.4n.filter <- subset(P1.4n, subset = Tnni3> 1)
P1.4n.filter <- subset(P1.4n, subset = Tnnt2>1)
P1.4n.filter <- subset(P1.4n, subset = Tnnc1>1)
P1.4n.filter <- subset(P1.4n, subset = Actc1> 1)
table(P1.4n.filter@meta.data$orig.ident)



#Merged Seurat objects.
Merged_seurat_object.P1 <- merge(x = P1.2n.filter, y = P1.4n.filter, add.cell.ids = c("P1.2n", "P1.4n"), project = "mouse_combined")
table(Merged_seurat_object.P1$orig.ident)



##-----------------------------------------------------Perform Integration----------------------------------------------------##
list <- SplitObject(Merged_seurat_object.P1, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)

#Perform the integration
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features, dims = 1:20)
combined.P1 <- IntegrateData(anchorset = anchors, dims = 1:20 )

combined.P1 <- ScaleData(combined.P1, vars.to.regress = "nCount_RNA", verbose = TRUE)

DefaultAssay(combined.P1) <- "integrated"

# Run the standard workflow for visualization and clustering

combined.P1 <- RunPCA(combined.P1, npcs = 30, verbose = FALSE)
combined.P1 <- RunUMAP(combined.P1, reduction = "pca", dims = 1:20)
combined.P1 <- FindNeighbors(combined.P1, reduction = "pca", dims = 1:20)
combined.P1 <- FindClusters(combined.P1, resolution = 0.15)




#-----------------------------------------------------------------------------------------------------------------------------------
#Read in a list of cell cycle markers

cc.genes <- readLines("")

#Aggregate the list into markers for G2/M and S phases
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
DefaultAssay(combined.P1) <- "RNA"

combined.P1 <- CellCycleScoring(combined.P1, s.features = s.genes, g2m.features = g2m.genes, 
                                set.ident = TRUE)
table(combined.P1$Phase)

##Alter threshold value
#combined.P1$Phase[combined.P1$S.Score < 0.1 & combined.P1$G2M.Score < 0.1] <- "G1"

table(combined.P1$Phase)


Idents(combined.P1) <- "Phase"


combined.P1@meta.data[["Phase"]] <- factor(combined.P1@meta.data[["Phase"]], levels = c("G1", "S", "G2M"))


p1 <- UMAPPlot(combined.P1,  pt.size = 1)


p2 <- DimPlot(combined.P1, reduction = "umap", split.by = "orig.ident", group.by = "Phase")



prop.table(table(Idents(combined.P1)))

table(Idents(combined.P1), combined.P1$orig.ident)

table(combined.P1$orig.ident)

d <- prop.table(table(Idents(combined.P1), combined.P1$orig.ident), margin = 2)

p3 <- tableGrob(format(round(d, 3)))

grid.arrange(p3)




tiff(file = "_CC_P1_CorrectSettings.tiff", width = 12000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1, p2, p3, ncol=3 ,widths=c(6,8,2)))

dev.off()
##############################################################



Idents(combined.P1) <- "seurat_clusters"

####################ClusterProfiler###########################

combined_markers.P1 <- FindAllMarkers(combined.P1, min.pct = 0.25, log2fc.threshold = 0.25)
head(combined_markers.P1)

write.csv2(combined_markers.P1, file="all_markers_P1.csv")


Num_Clusters <- unique(combined_markers.P1$cluster)

Cluster_gene_list <- list()
for (i in Num_Clusters){ 
  Cluster_genes <- subset(combined_markers.P1, combined_markers.P1$cluster == i & combined_markers.P1$avg_log2FC>0.2)
  Cluster_genes <- Cluster_genes %>% pull(gene)
  Cluster_gene_list[[i]] <- Cluster_genes
}




#### GO Analysis with Clusterprofiler
background <- row.names(combined.P1)

GO_list <- list()
GO_list_simplified <- list()
for (i in Num_Clusters){
  GO_Cluster <- enrichGO(gene  = Cluster_gene_list[[i]],
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         universe = background,
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.01)
  GO_list[[i]] <- as.data.frame(GO_Cluster)
  GO_list_simplified[[i]] <- as.data.frame(simplify(GO_Cluster, cutoff=0.6, by="p.adjust", select_fun=min))
  
}


top_GO <- list()

for (i in Num_Clusters){
  
  top_GO <- lapply(GO_list[["i"]], "[[",  "Description[[1]]")
  
}


cat(capture.output(print(GO_list), file="GO_list_P1.txt"))


top_GO <- list()

top_GO <- lapply(GO_list, `[[`, 2)

top_GO_first <- list()
top_GO_first <- lapply(top_GO, `[`, 1)


new.cluster.ids <- sapply(1:length(top_GO_first), function(i) top_GO_first[[i]])



new.cluster.ids <- make.names(new.cluster.ids,unique=T)


names(new.cluster.ids) <- levels(combined.P1)

combined.P1 <- RenameIdents(combined.P1, new.cluster.ids)


p1 <- UMAPPlot(combined.P1,  pt.size = 1)


p2 <- DimPlot(combined.P1, reduction = "umap", split.by = "orig.ident", group.by = "seurat_clusters")





prop.table(table(Idents(combined.P1)))

table(Idents(combined.P1), combined.P1$orig.ident)

table(combined.P1$orig.ident)

d <- prop.table(table(Idents(combined.P1), combined.P1$orig.ident), margin = 2)

p3 <- tableGrob(format(round(d, 2)))

grid.arrange(p3)



tiff(file = "_GOterms_P1.tiff", width = 13000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1 + NoLegend(), p2, p3, ncol=3 ,widths=c(6,10,6)))

dev.off()

