

library(Seurat)

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
combined <- merge(E165.2n, y = c(E165.4n, P1.2n, P1.4n, P5.2n, P5.4n), add.cell.ids = c("E165.2n", "E165.4n", "P1.2n", "P1.4n", "P5.2n", "P5.4n"), project = "mouse_combined")
table(combined$orig.ident)





#Prepare data
# Normalization
combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
# Most variable genes
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# scales data of all genes.
all.genes <-row.names(combined)
combined <- ScaleData(combined, features = all.genes, vars.to.regress = "nCount_RNA")

#Define markers
cardio <- c("Actc1", "Myh6", "Atp2a2", "Tnnc1", "Tnni3",  "Tnnt2")
fibro <- c("Col1a2", "Col3a1", "Vim", "Fstl1", "Gsn", "Fbln2", "Mmp2", "Ddr2", "Thy1", "Pdgfrb")
endo <- c( "Tie1", "Egfl7", "Pecam1", "Flt1", "Emcn", "Ednrb")
macro <- c( "Lgals3", "Cd68",  "Ptprc")
smo <- c( "Tagln", "Cald1", "Myh11", "Cnn1")


#Figure 1d:
tiff(file = "d_heatmap.tiff", width = 6500, height = 4000, units = "px", res = 600, pointsize = 12)

p1 <- DoHeatmap(combined, group.by = "orig.ident", features = c(cardio, fibro, endo, macro, smo), 
                slot = "data", disp.max = 5) + 
  scale_fill_gradientn(colors = c("#A6CEE3", "#1F78B4"))


p1 + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()










combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.15)


Idents(combined) <- "orig.ident"
FeaturePlot(object = combined,  split.by="orig.ident", features = c("E2f7", "E2f8"), order = T)
DotPlot(combined, features = c("E2f7", "E2f8"))


