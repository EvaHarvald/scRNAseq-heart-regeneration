
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



setwd("")

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






cc.genes <- readLines("")


#Aggregate the list into markers for G2/M and S phases
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
DefaultAssay(combined) <- "RNA"

combined.cc <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, 
                                set.ident = TRUE)
table(combined.cc$Phase)
##Alter threshold value
#combined.cc$Phase[combined.cc$S.Score < 0.1 & combined.cc$G2M.Score < 0.1] <- "G1"




Idents(combined.cc) <- "Phase"

# Define an order of cluster identities
my_levels <- c("G1", "G2M", "S")

# Relevel object@ident
combined.cc@meta.data[["Phase"]] <- factor(x = combined.cc@meta.data[["Phase"]], levels = my_levels)

p1 <- DimPlot(combined.cc, reduction = "umap",  group.by = "Phase")


p2 <- DimPlot(combined.cc, reduction = "umap", split.by = "orig.ident", group.by = "Phase")



prop.table(table(Idents(combined.cc)))

table(Idents(combined.cc), combined.cc$orig.ident)

table(combined.cc$orig.ident)

Idents(combined.cc) <- "Phase"

d <- prop.table(table(Idents(combined.cc), combined.cc$orig.ident), margin = 2)

p3 <- tableGrob(format(round(d, 3)))

grid.arrange(p3)



tiff(file = "Figure3a.tiff", width = 18000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1, p2, p3, ncol=3 ,widths=c(10,24,8)))

dev.off()





tiff(file = "_CellCycle_woTable.tiff", width = 20000, height = 4000, units = "px", res = 650, pointsize = 12)

grid.arrange(
  arrangeGrob(
    p1, p2, ncol=2 ,widths=c(6,24)))

dev.off()


Idents(combined.cc) <- "Phase"


tiff(file = "_CellCycle_dotplot.tiff", width = 4000, height = 2500, units = "px", res = 650, pointsize = 12)

plot <- DotPlot(combined.cc, cols = c("lightgrey", "darkblue"),
        features = c("Zeb1", "Aurka", "Ccne2","Ccng2", "Tead1", "Twist1","Ccnd1", "Cdkn1a",  "Ctgf"  ), dot.scale = 10)
plot + theme(axis.text.x = element_text(angle = 45,  hjust=1))

dev.off()


Idents(combined.cc) <- "orig.ident"



DotPlot(combined.cc, cols = c("lightgrey", "darkblue"),
        features = c("Zeb1",  "Ccne2", "Cdkn1a" ), dot.scale = 10)



Idents(combined.cc) <- "Phase"



DotPlot(combined.cc, cols = c("lightgrey", "darkblue"),
        features = c("Zeb1",  "Ccne2", "Cdkn1a" ), dot.scale = 10)


FeaturePlot(combined.cc, features = "Ccne2", split.by = "orig.ident", order = T)


FeaturePlot(combined.cc, features = "Ccne2", split.by = "Phase", order = T)

FeaturePlot(combined.cc, features = c("Ccne2", "Ccnd1"), split.by = "Phase", blend = T, order = T)

FeaturePlot(combined.cc, features = c("Ccne2", "Cdkn1a"), split.by = "Phase", blend = T, order = T)


FeaturePlot(combined.cc, features = c("Cdkn1a", "Zeb1"), split.by = "Phase", blend = T, order = T)

FeaturePlot(combined.cc, features = c("Ccnd1", "Zeb1"), split.by = "Phase", blend = T, order = T)

FeaturePlot(combined.cc, features = c("Ccne2", "Zeb1"), split.by = "Phase", blend = T, order = T)




DotPlot(combined.cc, features = c("Zeb1", "Ccng2",  "Ccne2", "Tead1", "Twist1", "Aurka", "Ccnd1", "Ctgf", "Cdkn1a", "Mstn"), dot.scale = 10)
