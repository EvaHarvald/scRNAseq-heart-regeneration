##### Figure1 + table

setwd("")
#E16.5 4n
dat <- Read10X(data.dir = "C:/Users/ebhar/Desktop/REPHEARTDATA/Mouse_old/E165_4n_old")
colnames(dat) = paste0("E16.5.4n_", colnames(dat) )
E16.5.4n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")

#Table in 1F:
Actc1 <- subset(E16.5.4n, subset = Actc1 > 1)
Tnni3 <- subset(E16.5.4n, subset = Tnni3 > 1)
Tnnc1 <- subset(E16.5.4n, subset = Tnnc1 > 1)
Tnnt2 <- subset(E16.5.4n, subset = Tnnt2 > 1)

#For table 1,e:
dim(Tnni3)
dim(Actc1)
dim(Tnnc1)
dim(Tnnt2)

#Preparing data and running PCA and tSNE:
E16.5.4n <- NormalizeData(object = E16.5.4n, normalization.method = "LogNormalize", 
                                        scale.factor = 10000)
# Most variable genes
E16.5.4n <- FindVariableFeatures(E16.5.4n, selection.method = "vst", nfeatures = 2000)



E16.5.4n <- ScaleData(E16.5.4n, vars.to.regress = "nCount_RNA", verbose = TRUE)



# Run the standard workflow for visualization and clustering

E16.5.4n <- RunPCA(E16.5.4n, npcs = 30, verbose = FALSE)
E16.5.4n <- RunUMAP(E16.5.4n, reduction = "pca", dims = 1:20)




tiff(file = "Tnni3_E1654n.tiff", width = 4500, height = 4000, units = "px", res = 600, pointsize = 12)
FeaturePlot(E16.5.4n, features = c("Tnni3"), 
            cols = c("lightgrey", "#1F78B4"), 
            ncol = 1, order =T, pt.size = 1)
dev.off()





rm(list=ls())


######################################################################


#E16.5 2n
dat <- Read10X(data.dir = "C:/Users/ebhar/Desktop/REPHEARTDATA/Mouse_old/E165_2n_old/mm10")
colnames(dat) = paste0("E16.5.2n_", colnames(dat) )
E16.5.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")


#Table in 1F:
Actc1 <- subset(E16.5.2n, subset = Actc1 > 1)
Tnni3 <- subset(E16.5.2n, subset = Tnni3 > 1)
Tnnc1 <- subset(E16.5.2n, subset = Tnnc1 > 1)
Tnnt2 <- subset(E16.5.2n, subset = Tnnt2 > 1)
#For table 1,e:
dim(Actc1)
dim(Tnni3)
dim(Tnnc1)
dim(Tnnt2)


#Preparing data and running PCA and tSNE:
E16.5.2n <- NormalizeData(object = E16.5.2n, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
# Most variable genes
E16.5.2n <- FindVariableFeatures(E16.5.2n, selection.method = "vst", nfeatures = 2000)



E16.5.2n <- ScaleData(E16.5.2n, vars.to.regress = "nCount_RNA", verbose = TRUE)



# Run the standard workflow for visualization and clustering

E16.5.2n <- RunPCA(E16.5.2n, npcs = 30, verbose = FALSE)
E16.5.2n <- RunUMAP(E16.5.2n, reduction = "pca", dims = 1:20)






tiff(file = "Tnni3_E1652n.tiff", width = 4500, height = 4000, units = "px", res = 600, pointsize = 12)
FeaturePlot(E16.5.2n, features = c("Tnni3"), 
            cols = c("lightgrey", "#1F78B4"), 
            ncol = 1, pt.size  = 1, order = T)
dev.off()






rm(list=ls())



################################


#P1 4n
dat <- Read10X(data.dir = "C:/Users/ebhar/Desktop/REPHEARTDATA/Mouse_old/P1_4n_old")
colnames(dat) = paste0("P1.4n_", colnames(dat) )
P1.4n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")


#Table in 1F:
Actc1 <- subset(P1.4n, subset = Actc1 > 1)
Tnni3 <- subset(P1.4n, subset = Tnni3 > 1)
Tnnc1 <- subset(P1.4n, subset = Tnnc1 > 1)
Tnnt2 <- subset(P1.4n, subset = Tnnt2 > 1)
#For table 1,e:
dim(Actc1)
dim(Tnni3)
dim(Tnnc1)
dim(Tnnt2)


#Preparing data and running PCA and tSNE:
P1.4n <- NormalizeData(object = P1.4n, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
# Most variable genes
P1.4n <- FindVariableFeatures(P1.4n, selection.method = "vst", nfeatures = 2000)


P1.4n <- ScaleData(P1.4n, vars.to.regress = "nCount_RNA", verbose = TRUE)



# Run the standard workflow for visualization and clustering

P1.4n <- RunPCA(P1.4n, npcs = 30, verbose = FALSE)
P1.4n <- RunUMAP(P1.4n, reduction = "pca", dims = 1:20)



tiff(file = "Tnni3_P1_4n.tiff", width = 4500, height = 4000, units = "px", res = 600, pointsize = 12)
FeaturePlot(P1.4n, features = c("Tnni3"), 
            cols = c("lightgrey", "#1F78B4"), 
            ncol = 1, order =T, pt.size = 1)
dev.off()






rm(list=ls())


######################################################################


#P1 2n
dat <- Read10X(data.dir = "C:/Users/ebhar/Desktop/REPHEARTDATA/Mouse_old/P1_2n_old")
colnames(dat) = paste0("P1.2n_", colnames(dat) )
P1.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")



#Table in 1F:
Actc1 <- subset(P1.2n, subset = Actc1 > 1)
Tnni3 <- subset(P1.2n, subset = Tnni3 > 1)
Tnnc1 <- subset(P1.2n, subset = Tnnc1 > 1)
Tnnt2 <- subset(P1.2n, subset = Tnnt2 > 1)
#For table 1,e:
dim(Actc1)
dim(Tnni3)
dim(Tnnc1)
dim(Tnnt2)


#Preparing data and running PCA and tSNE:
P1.2n <- NormalizeData(object = P1.2n, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
# Most variable genes
P1.2n <- FindVariableFeatures(P1.2n, selection.method = "vst", nfeatures = 2000)


P1.2n <- ScaleData(P1.2n, vars.to.regress = "nCount_RNA", verbose = TRUE)



# Run the standard workflow for visualization and clustering

P1.2n <- RunPCA(P1.2n, npcs = 30, verbose = FALSE)
P1.2n <- RunUMAP(P1.2n, reduction = "pca", dims = 1:20)




tiff(file = "Tnni3_P1_2n.tiff", width = 4500, height = 4000, units = "px", res = 600, pointsize = 12)
FeaturePlot(P1.2n, features = c("Tnni3"), 
            cols = c("lightgrey", "#1F78B4"), 
            ncol = 1, pt.size  = 1, order = T)
dev.off()





rm(list=ls())



###########################


#P5 4n
dat <- Read10X(data.dir = "C:/Users/ebhar/Desktop/REPHEARTDATA/Mouse_old/P5_4n_old")
colnames(dat) = paste0("P5.4n_", colnames(dat) )
P5.4n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")
#Table in 1F:
Actc1 <- subset(P5.4n, subset = Actc1 > 1)
Tnni3 <- subset(P5.4n, subset = Tnni3 > 1)
Tnnc1 <- subset(P5.4n, subset = Tnnc1 > 1)
Tnnt2 <- subset(P5.4n, subset = Tnnt2 > 1)

#For table 1,e:
dim(Actc1)
dim(Tnni3)
dim(Tnnc1)
dim(Tnnt2)


#Preparing data and running PCA and tSNE:
P5.4n <- NormalizeData(object = P5.4n, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
# Most variable genes
P5.4n <- FindVariableFeatures(P5.4n, selection.method = "vst", nfeatures = 2000)


P5.4n <- ScaleData(P5.4n, vars.to.regress = "nCount_RNA", verbose = TRUE)



# Run the standard workflow for visualization and clustering

P5.4n <- RunPCA(P5.4n, npcs = 30, verbose = FALSE)
P5.4n <- RunUMAP(P5.4n, reduction = "pca", dims = 1:20)



tiff(file = "Tnni3_P5_4n.tiff", width = 4500, height = 4000, units = "px", res = 600, pointsize = 12)
FeaturePlot(P5.4n, features = c("Tnni3"), 
            cols = c("lightgrey", "#1F78B4"), 
            ncol = 1, order =T, pt.size = 1)
dev.off()





rm(list=ls())


######################################################################


#P5 2n
dat <- Read10X(data.dir = "C:/Users/ebhar/Desktop/REPHEARTDATA/Mouse_old/P5_2n_old")
colnames(dat) = paste0("P5.2n_", colnames(dat) )
P5.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse-E16-5_2n")

#Table in 1F:
Actc1 <- subset(P5.2n, subset = Actc1 > 1)
Tnni3 <- subset(P5.2n, subset = Tnni3 > 1)
Tnnc1 <- subset(P5.2n, subset = Tnnc1 > 1)
Tnnt2 <- subset(P5.2n, subset = Tnnt2 > 1)

#For table 1,e:
dim(Actc1)
dim(Tnni3)
dim(Tnnc1)
dim(Tnnt2)


#Preparing data and running PCA and tSNE:
P5.2n <- NormalizeData(object = P5.2n, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
# Most variable genes
P5.2n <- FindVariableFeatures(P5.2n, selection.method = "vst", nfeatures = 2000)



P5.2n <- ScaleData(P5.2n, vars.to.regress = "nCount_RNA", verbose = TRUE)



# Run the standard workflow for visualization and clustering

P5.2n <- RunPCA(P5.2n, npcs = 30, verbose = FALSE)
P5.2n <- RunUMAP(P5.2n, reduction = "pca", dims = 1:20)



tiff(file = "Tnni3_P5_2n.tiff", width = 4500, height = 4000, units = "px", res = 600, pointsize = 12)
FeaturePlot(P5.2n, features = c("Tnni3"), 
            cols = c("lightgrey", "#1F78B4"), 
            ncol = 1, pt.size  = 1, order = T)
dev.off()





rm(list=ls())




