# Week 1 project work
# Installing packages (matrix)
library(Matrix)

# Alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')

# Loading single cell data for silique
a=load ("data/GSE226097_silique_230221.RData")

# Loading the original AraBOXcis network that was trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

# Examining the data
print(a)

dim(gbox)

length(clust)

# Looking at gbox variable
rownames(gbox)[1:10]

colnames(gbox)[1:10]

as.matrix(gbox[1:50, 1:3])

# Looking at clust variable
clust[1:8]

table(clust)

# Plotting a graph using the sort function
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Silique')

dim(araboxcis)

araboxcis[1:4,]

# Plotting a histogram 
hist(araboxcis[,3])

tfs=unique(araboxcis[,1])

tfSubs=tfs[which(tfs %in% rownames(gbox))]

length(tfSubs)

# Filtering the genes and cells with very low values
dim(gbox)

thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

# Visualizing with UMAP
# Installing packages (Umap)
install.packages("umap")
library(umap)

# Using UMAP to visualize
gbox.umap <- umap(gbox_filtered)

# Do all cell type clusters group together if only looking at g-box related genes
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Silique', xlab='UMAP Component 1', ylab='UMAP Component 2')

pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Silique', xlab='UMAP Component 1', ylab='UMAP Component 2')
