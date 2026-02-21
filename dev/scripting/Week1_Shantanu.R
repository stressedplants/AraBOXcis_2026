library(Matrix)
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

source('dev/utilities/dataprocessingHelperFunctions.R')

a=load('data/GSE226097_stem_230221.RData')
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

print(a)

# Getting the size of the matrix (gbox), number or rows and columns
dim(gbox)

# Getting the length of the vetcor(clust) 
length(clust)

# Accessing and printing the first 10 row names in gbox
# Rownames are names of genes 
rownames(gbox)[1:10]

# Printing only the 3rd row name in gbox
rownames(gbox)[3]

# Accessing and printing the 2nd and 4th row name only
rownames(gbox)[c(2, 4)] 

# Accessing the first 10 columns names 
# Column names are individual cells (unique IDs)
colnames(gbox)[1:10]

# Accessing the 4th column name
colnames(gbox)[4]

# Accessing the first 50 rows and first 3 columns of the gbox matrix
as.matrix(gbox[1:50, 1:3])

# accessing row 1 and all columns 
as.matrix(gbox[1,])

# accessing column 1 and all rows
as.matrix(gbox[,1])

# Printing cluster designation of first 8 cells
# Each cell is designated a cluster ranging from 1-15
clust[1:8]

# Getting a count of the number of cells in each cluster 
table(clust)

# Plotting a frequency graph of the clust table 
plot(table(sort(clust)), xlab = 'Cluster Name', ylab= ' Number of Cells', main='Stem')

# Checking the size of the araboxcis data (50000 rows, 3 columns)
# Column 1= Transcription factors, Column 2= Target gene and Column 3= Regulatory edge score
dim(araboxcis)

# Accessing the first 4 rows 
araboxcis[1:4,]

# Plotting histogram of column 3 (regulatory edge scores) in araboxis data
hist(araboxcis[,3])

# Re-plotting histogram of column 3 (regulatory edge scores) in araboxis data
hist(araboxcis$score, xlab = 'Regulatory Edge Score', ylab = 'Frequency', main = 'Distribution of Regulatory Edge Scores')

# Getting a list of the unique transcription factors in column 1
tfs=unique(araboxcis[,1])
# Checking how many unique transcription factors are present in araboxis dataset
length(tfs)

# Creating a subset/ filtering the gbox data to only include transcription factors 
# that are present in the araboxis dataset
tfSubs=tfs[which(tfs %in% rownames(gbox))]

# Checking how many common TF's are present between the gbox and araboxis datasets 
# (unique tf's- Stem:139)
length(tfSubs)

## which()- gives position of the value in data

dim(gbox)

# Getting rid of cells that have less than 1% of genes expressed in them

# Setting threshold to 0.01 or 1%
thresh=0.01
# For each column (Cell) iterating through the rows (genes) and counting how many rows(genes) have a value greater than 0
# For each column(cell), counting the number of rows(genes) with value a greater than 0 
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})

# In numberGenesPerCell 
# including columns (cells) that have a count (number of genes>0) greater than the threshold (0.01 x 2096=20.96)
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

# Creating a new matrix that keeps the same rows (genes) but only includes columns in includeCells
gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

# Method 2
numberGenesPerCell <- colSums(gbox > 0)
includeCells <- which(numberGenesPerCell > (thresh * nrow(gbox)))
gbox_filtered <- gbox[, includeCells]
dim(gbox_filtered)

# Getting rid of genes that are expressed in less than 1% of cells 
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)


# Method 2 
# numberCellsPerGene <- rowSums(gbox > 0)
# includeGenes <- which(numberCellsPerGene > (thresh * ncol(gbox)))
# gbox_filtered <- gbox_filtered[includeGenes,]
# dim(gbox_filtered)

install.packages('umap')
library(umap)
library(viridis)
#Let us visualise it using UMAP
gbox.umap <- umap(gbox_filtered)

#Do the cell type clusters group together if we only look at G-box related genes?
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], 
     pch=20, main='UMAP Stem', xlab='UMAP Component 1', ylab='UMAP Component 2')

legend("top",
       legend = unique(clust[includeCells]),
       fill = colours,
       title = "Clusters",
       horiz = TRUE,
       cex=0.45)

## Testing different n_neighbors for the UMAP

# Creating a 2x3 grid of plots
par(mfrow = c(2, 3)) 

# Creating a list of n_neighbors to test
neighbors <- c(2, 3, 4, 5, 6)

for (n in neighbors) {
  config <- umap.defaults
  config$n_neighbors <- n
  gbox_neighbors<- umap(gbox_filtered, config = config)
  
  plot(gbox_neighbors$layout[,1], gbox_neighbors$layout[,2], col=colours[clust[includeCells]], 
       pch=20,
       main=paste0("n_neighbors:", n))
}

## Testing different min_dist for the UMAP
par(mfrow = c(2, 3))

distance <- c(0.01,0.1, 0.2, 0.4, 0.6, 0.8)

for (d in distance) {
  config <- umap.defaults
  config$min_dist <- d
  gbox_dist<- umap(gbox_filtered, config = config)
  
  plot(gbox_dist$layout[,1], gbox_dist$layout[,2], col=colours[clust[includeCells]], 
       pch=20,
       main=paste0("min_dist:", d))
}

# Creating a pca plot
pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]],
     pch=20, main='PCA UMAP Stem', xlab='UMAP Component 1', ylab='UMAP Component 2')




