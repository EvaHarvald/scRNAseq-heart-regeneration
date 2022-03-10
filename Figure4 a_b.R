
library(DOSE)
library(enrichplot)
library(ggnewscale)


library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)
library(patchwork)
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
library(metap)
library(gridExtra)
library(raster)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(clusterProfiler)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")



setwd("T:/CVR/LMCC/REPHEART 2021/Figure 4/Files from EH")

#Load all six samples
dat <- Read10X(data.dir = "")
colnames(dat) = paste0("E16.5.2n_", colnames(dat) )
E165.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")

dat1 <- Read10X(data.dir = "")
colnames(dat1) = paste0("E16.5.4n_", colnames(dat1) )
E165.4n <- CreateSeuratObject(counts = dat1, min.cells = 5, min.features = 200, project = "E165_4n")

dat2 <- Read10X(data.dir = "")
colnames(dat2) = paste0("P1.2n_", colnames(dat2) )
P1.2n <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "P1_2n")

dat3 <- Read10X(data.dir = "")
colnames(dat3) = paste0("P1.4n_", colnames(dat3) )
P1.4n <- CreateSeuratObject(counts = dat3, min.cells = 5, min.features = 200, project = "P1_4n")

dat4 <- Read10X(data.dir = "")
colnames(dat4) = paste0("P5.2n_", colnames(dat4) )
P5.2n <- CreateSeuratObject(counts = dat4, min.cells = 5, min.features = 200, project = "P5_2n")

dat5 <- Read10X(data.dir = "")
colnames(dat5) = paste0("P5.4n_", colnames(dat5) )
P5.4n <- CreateSeuratObject(counts = dat5, min.cells = 5, min.features = 200, project = "P5_4n")








#Merged Seurat objects.
combined <- merge(E165.2n, y = c(E165.4n, P1.2n, P1.4n, P5.2n, P5.4n), add.cell.ids = c("E165.2n", "E165.4n", "P1.2n", 
                                                                                        "P1.4n", "P5.2n", "P5.4n"), project = "mouse_combined")
table(combined$orig.ident)


##-----------------------------------------------------Filtering cells based on markers----------------------------------------------------##
#remove non CMs
combined <- subset(combined, subset = Tnni3> 1)
combined <- subset(combined, subset = Tnnt2>1)
combined <- subset(combined, subset = Tnnc1>1)
combined <- subset(combined, subset = Actc1> 1)
table(combined@meta.data$orig.ident)



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")




# normalize and identify variable features for each dataset independently

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined, vars.to.regress = "nCount_RNA", verbose = TRUE)


# Run the standard workflow for visualization and clustering
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.15)

#############################################################



Idents(combined) <- "seurat_clusters"

####################ClusterProfiler###########################

combined_markers <- FindAllMarkers(combined, min.pct = 0.25, log2fc.threshold = 0.25)
head(combined_markers)

write.csv2(combined_markers, file="all_markers.csv")


Num_Clusters <- unique(combined_markers$cluster)

Cluster_gene_list <- list()
for (i in Num_Clusters){ 
  Cluster_genes <- subset(combined_markers, combined_markers$cluster == i & combined_markers$avg_log2FC>0.5)
  Cluster_genes <- Cluster_genes %>% pull(gene)
  Cluster_gene_list[[i]] <- Cluster_genes
}




#### GO Analysis with Clusterprofiler
background <- row.names(combined)

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



top_GO <- list()

top_GO <- lapply(GO_list, `[[`, 2)

top_GO_first <- list()
top_GO_first <- lapply(top_GO, `[`, 1)


new.cluster.ids <- sapply(1:length(top_GO_first), function(i) top_GO_first[[i]])



new.cluster.ids <- make.names(new.cluster.ids,unique=T)


names(new.cluster.ids) <- levels(combined)

combined <- RenameIdents(combined, new.cluster.ids)











##########################
#Read in a list of cell cycle markers

cc.genes <- readLines("")


#Aggregate the list into markers for G2/M and S phases
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
DefaultAssay(combined) <- "RNA"

combined.cc <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, 
                                set.ident = TRUE)
table(combined.cc$Phase)

#Figure a:
UMAP(combined.cc)



##########################
#Compare "cell dividion" clusters
Idents(combined.cc) <- "Phase"


subset.G2M <- subset(x = combined.cc, idents =  "G2M")




table(subset.G2M$orig.ident)


Idents(subset.G2M) <- "orig.ident"
subset.G2M.4n <- subset(x = subset.G2M, idents =  c("E16.5.4n", "P5.4n"))


subset.G2M.4n <- ScaleData(subset.G2M.4n, vars.to.regress = "nCount_RNA", verbose = TRUE)


# Run the standard workflow for visualization and clustering
subset.G2M.4n <- RunPCA(subset.G2M.4n, npcs = 30, verbose = FALSE)
subset.G2M.4n <- RunUMAP(subset.G2M.4n, reduction = "pca", dims = 1:20)
subset.G2M.4n <- FindNeighbors(subset.G2M.4n, reduction = "pca", dims = 1:20)
subset.G2M.4n <- FindClusters(subset.G2M.4n, resolution = 0.15)

Idents(subset.G2M.4n) <- "orig.ident"

#Figure b:

tiff(file = "UMAP_subset.tiff", width = 4000, height = 3000, units = "px", res = 650, pointsize = 12)
UMAPPlot(subset.G2M.4n)
dev.off()


table(subset.G2M.4n$orig.ident)



cluster.markers.subset.G2M <- FindMarkers(subset.G2M.4n, ident.1 = "E16.5.4n", ident.2 = "P5.4n", min.pct = 0.25)


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


#Figure d:

tiff("ccGOs_upincluster3_G2M_rainbow.tiff",width=6000,height=2000, units = "px", res = 350)
cnetplot(GO_Cluster,  circular = TRUE, colorEdge = TRUE,
         showCategory = c("mitotic nuclear division","cell division", "chromosome condensation", "regulation of mitotic cell cycle", "nuclear division"))
dev.off()


write.csv2(GO_list, file="GO_list_upInCluster3_celldivision_G2M4ncells.csv")








setDT(cluster.markers.subset.G2M, keep.rownames = "Genes")



top <- cluster.markers.subset.G2M  %>% top_n(n = 100, wt = avg_log2FC)


#Figure c:
tiff("heatmap.tiff",width=2000,height=3200, units = "px", res = 300)
DoHeatmap(subset.G2M.4n, group.by = "ident", features = top$Genes, 
          group.colors = c("lightgrey", "black"), slot = "scale.data") +
  theme(axis.text = element_text(size = 10)) +
  scale_fill_gradientn(colors = c("white", "aliceblue", "steelblue4"), na.value = "white") +
  theme(legend.position ="top",legend.direction = "horizontal")
dev.off()

write.csv2(cluster.markers.subset.G2M, file="cluster.markers_celldivision.subset.G2M.4n.csv")










dt <- read.delim("TFs_upincluster3_possum_G2M.txt", header=T)

dt <- dt[dt$Fisher.score > 0.5,]
dt$TF <- factor(dt$TF, levels = dt$TF[order(dt$Fisher.score)])


dt$TF <- factor(dt$TF, levels = dt$TF[(order(dt$Fisher.score))])

#Figure e:
tiff("TFs_upincluster3_possum_G2M.tiff",width=1400,height=1200, units = "px", res = 300)
ggplot(dt[1:20,], aes(y=TF, x=Fisher.score)) + 
  geom_point(stat="identity")+
  geom_point(aes(colour = Target.gene.hits)) +
  theme(legend.position="right")+
  guides(size = guide_legend(order = 1))+
  theme_bw()
dev.off()

##################################################


