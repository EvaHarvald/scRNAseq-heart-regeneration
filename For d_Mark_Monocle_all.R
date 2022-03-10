

library(Seurat)
library(monocle)
library(ggplot2)

packageVersion("monocle")
packageVersion("Seurat")
# Load samples - change directory

dat <- Read10X(data.dir = "")
colnames(dat) = paste0("E16.5.2n_", colnames(dat) )
E16.5.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse.E16.5_2n")




dat <- Read10X(data.dir = "")
colnames(dat) = paste0("P1.2n_", colnames(dat) )
P1.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse.P1_2n")




dat <-  Read10X(data.dir = "")
colnames(dat) = paste0("P5.2n_", colnames(dat) )
P5.2n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse.P5_2n")



dat <- Read10X(data.dir = "")
colnames(dat) = paste0("E16.5.4n_", colnames(dat) )
E16.5.4n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse.E16.5_4n")




dat <- Read10X(data.dir = "")
colnames(dat) = paste0("P1.4n_", colnames(dat) )
P1.4n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse.P1_4n")




dat <- Read10X(data.dir = "")
colnames(dat) = paste0("P5.4n_", colnames(dat) )
P5.4n  <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "mouse.P5_4n")







#Combine: 

combined.all <- merge(x = E16.5.2n, y = c(P1.2n, P5.2n, E16.5.4n, P1.4n, P5.4n), add.cell.ids = c("E16.5.2n","P1.2n", "P5.2n", "E16.5.4n", "P1.4n", "P5.4n"), 
                       project = "mouse2n")


table(combined.all@meta.data$orig.ident)

combined.all <- subset(combined.all, subset = Tnni3>1)
combined.all <- subset(combined.all, subset = Tnnt2>1)
combined.all <- subset(combined.all, subset = Tnnc1>1)
combined.all <- subset(combined.all, subset = Actc1> 1)
table(combined.all@meta.data$orig.ident)

#Assign cell cycle phase - Load regev_lab_cell_cycle_genes.txt

cc.genes <- readLines(con = "")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

combined.all <- CellCycleScoring(object = combined.all, s.features = s.genes, g2m.features = g2m.genes, 
                                 set.ident = TRUE)



#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(combined.all@assays$RNA@data), 'sparseMatrix')

pd <- new("AnnotatedDataFrame", data = combined.all@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
combined.monocle <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
head(pData(combined.monocle))


class(combined.monocle)
slotNames(combined.monocle)

#This function estimates the size factors using the "median ratio method"
combined.monocle <- estimateSizeFactors(combined.monocle)

#This function obtains dispersion estimates for a count data set. For each condition (or collectively for all 
#conditions, see 'method' argument below) it first computes for each gene an empirical dispersion value (a.k.a. a 
#raw SCV value), then fits by regression a dispersion-mean relationship and finally chooses for each gene a dispersion 
#parameter that will be used in subsequent tests from the empirical and the fitted value according to the 'sharingMode' argument.
combined.monocle <- estimateDispersions(combined.monocle, method = "pooled")

combined.monocle

#Sets the global expression detection threshold to be used with this CellDataSet. 
#Counts how many cells each feature in a CellDataSet object that are detectably expressed above a minimum threshold. 
#Also counts the number of genes above this threshold are detectable in each cell.
combined.monocle <- detectGenes(combined.monocle, min_expr = 0.1)
print(head(fData(combined.monocle)))


head(pData(combined.monocle))
summary(pData(combined.monocle)$num_genes_expressed)
summary(fData(combined.monocle)$num_cells_expressed)


# standardise to Z-distribution
x <- pData(combined.monocle)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
summary(x_1)


df <- data.frame(x = x_1)
ggplot(df, aes(x)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')




ggplot(pData(combined.monocle), aes(num_genes_expressed, nCount_RNA, color = orig.ident)) + geom_point()
ggplot(pData(combined.monocle), aes(num_genes_expressed, nCount_RNA, color = Phase)) + geom_point()


#Removing cells with too low or too high num_genes_expressed ##

valid_cells <- row.names(subset(pData(combined.monocle),
                                num_genes_expressed < 3000 & num_genes_expressed > 500))

combined.monocle <- combined.monocle[,valid_cells]



#The first step is determining a subset of genes to use for clustering; this is because not all genes are informative,
#such as those that are lowly expressed. The approach is to select gene based on their average expression and variability across cells.
#The dispersionTable() function calculates the mean and dispersion values.
# Retrieve table of values specifying the mean-variance of genes
disp_table <- dispersionTable(combined.monocle)
head(disp_table)

#We will select genes, which have a mean expression >= 0.1, to use in the clustering step. 
#The setOrderingFilter() function allows us to indicate which genes we want to use for clustering. 
#The plot_ordering_genes() function plots mean expression against the empirical dispersion and highlights the set of genes (as black dots) that will be used for clustering.

table(disp_table$mean_expression>=0.1)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

#The function marks genes that will be used for clustering in subsequent calls to clusterCells. The list of selected genes can be altered at any time.
combined.monocle <- setOrderingFilter(combined.monocle, unsup_clustering_genes$gene_id)
plot_ordering_genes(combined.monocle)

head(fData(combined.monocle))

plot_pc_variance_explained(combined.monocle, return_all = FALSE)


#############################################

# Cluster for tSNE
combined.monocle <- reduceDimension(combined.monocle, max_components = 2, num_dim = 8,
                                    reduction_method = 'tSNE', verbose = TRUE)
# perform unsupervised clustering requesting 15-1 clusters
combined.monocle <- clusterCells(combined.monocle, num_clusters = 15)
head(pData(combined.monocle))


tiff(file = "tSNE_combined.tiff", width = 5500, height = 7500, units = "px", res = 950, pointsize = 12)
plot_cell_clusters(combined.monocle)
dev.off()

tiff(file = "tSNE_combined_orig_ident.tiff", width = 5500, height = 6500, units = "px", res = 950, pointsize = 12)
plot_cell_clusters(combined.monocle, color_by = 'orig.ident')
dev.off()

tiff(file = "tSNE_combined_CCphase.tiff", width = 5500, height = 6500, units = "px", res = 950, pointsize = 12)
plot_cell_clusters(combined.monocle, color_by = 'Phase')
dev.off()

tiff(file = "tSNE_combined_myomarkers.tiff", width = 11000, height = 6500, units = "px", res = 950, pointsize = 12)
plot_cell_clusters(combined.monocle, color_by = 'orig.ident', markers = c("Myh6", "Myh7"))
dev.off()
#Myh6 increases over time, whereas Myh7 decreases

###################################################

# Choose genes
expressed_genes <- row.names(subset(fData(combined.monocle), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(combined.monocle[expressed_genes,],
                                      fullModelFormulaStr = "~orig.ident")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.0001))


combined.monocle <- setOrderingFilter(combined.monocle, ordering_genes)
plot_ordering_genes(combined.monocle)



#Reduce dimensions

combined.monocle <- reduceDimension(combined.monocle, max_components = 2,
                                    method = 'DDRTree')

#Order cells in the trajectory returning Pseudotime and State
#Add roote_state here according to E16.5 being the start
combined.monocle <- orderCells(combined.monocle)


#Set wd 

tiff(file = "trajectory_all.TFs.tiff", width = 5500, height = 4500, units = "px", res = 950, pointsize = 12)
plot_cell_trajectory(combined.monocle, markers = c("Arnt", "Myc",
                                                      "Mycn","Zeb1", "Sp1"), use_color_gradient = TRUE)
dev.off()








GM_state <- function(combined.monocle){
  if (length(unique(pData(combined.monocle)$State)) > 1){
    T0_counts <- table(pData(combined.monocle)$State, pData(combined.monocle)$orig.ident)[, "E16.5.4n"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}


combined.monocle <- orderCells(combined.monocle, root_state = GM_state(combined.monocle))
plot_cell_trajectory(combined.monocle, color_by = "Pseudotime")





tiff(file = "trajectory_all.orig.ident.tiff", width = 5500, height = 4500, units = "px", res = 950, pointsize = 12)
plot_cell_trajectory(combined.monocle, color_by = "orig.ident")
dev.off()

#"State" is just Monocle's term for the segment of the tree

tiff(file = "trajectory_all.state.tiff", width = 5500, height = 4500, units = "px", res = 950, pointsize = 12)
plot_cell_trajectory(combined.monocle, color_by = "State")
dev.off()

tiff(file = "trajectory_all.pseudotime.tiff", width = 5500, height = 4500, units = "px", res = 950, pointsize = 12)
plot_cell_trajectory(combined.monocle, color_by = "Pseudotime")
dev.off()

tiff(file = "trajectory_all.phase.tiff", width = 5500, height = 4500, units = "px", res = 950, pointsize = 12)
plot_cell_trajectory(combined.monocle, color_by = "Phase")
dev.off()




#Finding Genes that Change as a Function of Pseudotime

my_pseudotime_de <- differentialGeneTest(combined.monocle,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 8)


library(dplyr)
my_pseudotime_de %>% arrange(qval) %>% head()

# save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$gene_short_name

plot_genes_in_pseudotime(combined.monocle[my_pseudotime_gene,])



# cluster the top 10 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(10) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$gene_short_name

my_pseudotime_cluster <- plot_pseudotime_heatmap(combined.monocle[gene_to_cluster,],
                                                 num_clusters = 3,
                                                 cores = 8,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)



my_cluster <- cutree(my_pseudotime_cluster$tree_row, 3)
my_cluster
# genes in cluster 1
my_pseudotime_de[names(my_cluster[my_cluster == 1]),"gene_short_name"]

#The BEAM() function takes a CellDataSet that has been ordered with orderCells() and a branch point in the trajectory.
#A table of genes is returned with significance values that indicate whether genes have expression patterns that are 
#branch dependent.



res <- length(combined.monocle@auxOrderingData[[combined.monocle@dim_reduce_type]]$branch_points)

for (i in 1:res) {
  BEAM_res <- BEAM(combined.monocle, branch_point = i, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval", "num_cells_expressed")]
  
  my_branched_heatmap <- plot_genes_branched_heatmap(combined.monocle[row.names(subset(BEAM_res, qval < 1e-4)),],
                                                     branch_point = i,
                                                     num_clusters = 4,
                                                     cores = 1,
                                                     use_gene_short_name = TRUE,
                                                     show_rownames = TRUE,
                                                     return_heatmap = TRUE)
  
  filename <- paste("Branch point", i, ".txt", sep="")
  write.csv2(BEAM_res, filename)
  
}

plot_list = list()
for (i in 1:res) {
  p = plot_genes_branched_heatmap(combined.monocle[row.names(subset(BEAM_res, qval < 1e-4)),],
                                                     branch_point = i,
                                                     num_clusters = 4,
                                                     cores = 1,
                                                     use_gene_short_name = TRUE,
                                                     show_rownames = TRUE,
                                                     return_heatmap = TRUE)
  
  plot_list[[i]] = p
  
}



for (i in 1:res) {
  file_name = paste("heatmap_BP", i, ".tiff", sep="")
  tiff(file_name, width = 5500, height = 7500, units = "px", res = 950, pointsize = 12)
  print(plot_list[[i]])
  dev.off()
}








