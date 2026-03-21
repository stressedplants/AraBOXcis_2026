#week 1
install.packages('Matrix')
library(Matrix)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

source('dev/utilities/dataprocessingHelperFunctions.R')

a=load('data/GSE226097_rosette_21d_230221.RData')
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)
print(a)

dim(gbox)
length(clust)
rownames(gbox)[1:10]
colnames(gbox)[1:10]
as.matrix(gbox[1:50, 1:3])
clust[1:8]
table(clust)

plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Rosette_21d')

dim(araboxcis)
araboxcis[1:4,]

hist(araboxcis[,3])
tfs=unique(araboxcis[,1])
tfSubs=tfs[which(tfs %in% rownames(gbox))]
length(tfSubs)
#Rosette_21d:145

thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

install.packages('umap')
library(umap)
gbox.umap <- umap(gbox_filtered)
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Rosette_21d', xlab='UMAP Component 1', ylab='UMAP Component 2')

pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Rosette_21d', xlab='UMAP Component 1', ylab='UMAP Component 2')

#week 3
install.packages(c("Matrix","Rcpp","rlang","lattice","openssl"), force = TRUE)
BiocManager::install("GENIE3", force = TRUE)
library(Matrix)
library(GENIE3)

net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=5)
save(net, file='Rosette_21d_network_nTree_5.RData')
ginieOutput=convertToAdjacency(net, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]

a=load('Rosette_21d_network_nTree_5.RData')
newNet=GENIE3::getLinkList(net)
top1000 = newNet[1:1000, ]
colnames(top1000) = c("TF", "target", "score")
top1000$dataset = "Rosette_21d"
geneCounts = table(c(top1000$TF, top1000$target))
validGenes = names(geneCounts[geneCounts > 1])
top1000_filtered = top1000[
  top1000$TF %in% validGenes &
    top1000$target %in% validGenes, ]
head(top1000_filtered)
dim(top1000_filtered)
View(top1000)
write.csv(top1000, "top1000_edges.csv", row.names = FALSE)

araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

genesInNet=unique(c(newNet[,1], newNet[,2]))
araboxcisFiltered=araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]
newNetTopEdges=newNet[1:length(araboxcisFiltered[,1]),]
edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')
length(which(edgesNew %in% edgesOld))
length(which(! (edgesNew %in% edgesOld)))
length(which(! (edgesOld %in% edgesNew)))

tfsNew=table(newNetTopEdges[,1])
tfsOld=table(araboxcisFiltered[,1])[names(tfsNew)]
hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs')
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')
sort(tfsNew, decreasing=TRUE)[1:20]

install.packages('igraph')
library(igraph)
install.packages('network')
library(network)
install.packages('pheatmap')
library(pheatmap)

simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))
node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:20]
plot(sort(node_betweenness))

node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]
plot(sort(node_hub))
plot(node_hub_all, node_centrality_all)
plot(node_hub_all, node_betweenness_all)

a=load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')
pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)

atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)
atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)
pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]
pheatmap(pafwayInterestingOnly)
pheatmap(log(pafwayInterestingOnly, 10))
