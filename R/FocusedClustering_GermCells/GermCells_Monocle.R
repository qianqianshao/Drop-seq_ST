### R script for pseudotemporal ordering of germ cells using Monocle in Nov 2017 by Qianyi
### Related to Figure S2C right panel: germ cells (N=20,646) with >1k detected genes 

### Read Data
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
file="figJul2017_10CellTypes_MouseAdultST24mergedclusters/monocle_"
load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))
dge20=dge   # Gene Filter 1 (N=24482): genes >20totalUMIsOverAllCells & >15Cells +31 genes
table(dge@data.info[,52:53]) 
# Cluster31SeriationOLO_GeneFilter1 is the reordered cluster ID 1-31
# Cluster 8-31 are germ cells
### Extract germ cells
cells.use=rownames(dge@data.info)[which(dge@data.info$Cluster31SeriationOLO_GeneFilter1>=8)]
anno <- data.frame(V1=rownames(dge@data.info),V2=dge@data.info$Cluster31SeriationOLO_GeneFilter1,dge@data.info)[which(dge@data.info$Cluster31SeriationOLO_GeneFilter1>=8),]
table(anno$V2)
summary(anno$nGene)
summary(anno$nUMI)
anno=anno[,-c(10:48,65:66)]
for(i in 1:24){
    anno$V2[which(anno$Cluster31SeriationOLO_GeneFilter1 == i+7)] <- i
}
anno$V2=factor(anno$V2,levels=1:24)
table(anno$V2)
### Remove germ cells with <=1000 Genes, because germ cells tend to have more transcripts compared with somatic cells
anno=anno[anno$nGene > 1000,]
table(anno$V2)
dim(anno)  # 20646 cells
all.col=col.use[anno$V2]

## ----package_loads, include=FALSE, eval=TRUE-----------------------------
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)

knitr::opts_chunk$set(autodep=TRUE, cache=FALSE, warning=FALSE)
set.seed(0)

## ----init_monocle, include=FALSE, eval=TRUE------------------------------
library(HSMMSingleCell)
library(monocle)
data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)

###### 1. Store data in CellDataSet
## ----build_cell_data_Set_RPC, eval=TRUE----------------------------------
umi_matrix=dge@raw.data[,anno$V1]
sample_sheet=anno
gene_annotation=data.frame(gene_short_name=rownames(umi_matrix),num_cells_expressed=rowSums(as.matrix(umi_matrix)>0))
rownames(gene_annotation)=rownames(dge@data)

dim(umi_matrix) # [1] 24482 20646
dim(sample_sheet) # 20646  25
dim(gene_annotation)
summary(gene_annotation[,2]) # min=0; max=20240

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

# First create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(as(umi_matrix, "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,,
                       expressionFamily=negbinomial.size())

## ----estimate_size_and_dispersion, eval=TRUE-----------------------------
HSMM <- estimateSizeFactors(HSMM)
# note: change add_rownames to tibble::rownames_to_column() in estimateDispersionsForCellDataSet function
HSMM <- estimateDispersions(HSMM)
## Removing 185 outliers

object=HSMM
modelFormulaStr="~ 1"; relative_expr=TRUE; min_cells_detected=1; remove_outliers=TRUE; cores=1


###### 2. Filtering low-quality cells (recommended)
## ----detect_genes, eval=TRUE---------------------------------------------
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
# keep genes expressed in >=10 cells
length(expressed_genes) # 24140

## ----show_pData, eval = TRUE---------------------------------------------
print(head(pData(HSMM)))

## ----select_cells, eval = FALSE------------------------------------------
#  valid_cells <- row.names(subset(pData(HSMM), Cells.in.Well == 1 & Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1000000))
#  HSMM <- HSMM[,valid_cells]

### double-check the expression follows a roughly-normal distribution
## ----lognormal_plot, eval=TRUE, fig.width = 3, fig.height = 2, fig.align="center"----
# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
pdf(paste0(file,"distribution_StdGeneExp.pdf"))
qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') +
  xlab("Standardized log(StdGeneExp)") +
  ylab("Density")
dev.off()


######### 3. Constructing Pseudotime Trajectories
#### step 1: Choose genes that define progress
## select genes that differ between clusters (recommended) 
## clustering cells using genes selected by unsupervised "dpFeature"
## select superset of feature genes as genes expressed in at least 5% of all cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
fData(HSMM)$use_for_ordering <- fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)
table(fData(HSMM)$use_for_ordering)
#FALSE  TRUE 
#14685  9797
HSMM <- setOrderingFilter(HSMM, fData(HSMM)[which(fData(HSMM)$use_for_ordering),]$gene_short_name)
pdf(paste0(file,"2.choosegenes_forclustering_24GermClusters.pdf"))
plot_ordering_genes(HSMM)
dev.off()
### perform PCA
# change default nv=pcs.ompute=20, maxit=100, fastpath=TRUE, work=nv+7 to maxit=500,work=500,fastpath=FALSE for irlba
pdf(paste0(file,"2.PCAvariance_24GermClusters.pdf"))
plot_pc_variance_explained(HSMM, return_all = F,norm_method = "log")
dev.off()

### perform tSNE
HSMM <- reduceDimension(HSMM, max_components = 2, norm_method = 'log',num_dim = 5, reduction_method = 'tSNE', verbose = T)
#Remove noise by PCA ..., reduction_method = 'tSNE', verbose = T)
#Reduce dimension by tSNE ...

### Cluster Cells
HSMM <- clusterCells(HSMM, verbose = F,num_clusters=25) 
# num_clusters = Number of clusters + 1
#Distance cutoff calculated to 6.781424 
#Warning messages:
#1: In if (method == "DDRTree") { :
#  the condition has length > 1 and only the first element will be used
#2: In if (method == "densityPeak") { :
#  the condition has length > 1 and only the first element will be used

table(pData(HSMM)$V2,pData(HSMM)$Cluster)
pdf(paste0(file,"2.ReCluster_24GermClusters.pdf"),width=5,height=6.5)
plot_cell_clusters(HSMM, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM, color_by = 'as.factor(V2)')
dev.off()

plot_rho_delta(HSMM, rho_threshold = 2, delta_threshold = 4 )
HSMM <- clusterCells(HSMM,
                         rho_threshold = 2,
                         delta_threshold = 4,
                         skip_rho_sigma = T,
                         verbose = F)
plot_cell_clusters(HSMM, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM, color_by = 'as.factor(Hours)')

### select genes that differ between clusters
HSMM_expressed_genes <-  row.names(subset(fData(HSMM),num_cells_expressed >= 10))
length(HSMM_expressed_genes)
#[1] 24140

clustering_DEG_genes <- differentialGeneTest(HSMM[HSMM_expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 10)                              
HSMM_ordering_genes <-
      row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
HSMM <- setOrderingFilter(HSMM, ordering_genes = HSMM_ordering_genes)
length(HSMM_ordering_genes) # 1000

pdf(paste0(file,"2.choosegenes_usedmarkers_24GermClusters.pdf"))
plot_ordering_genes(HSMM)
dev.off()

####step 2: reduce data dimensionality
HSMM <- reduceDimension(HSMM, method = 'DDRTree')

####step 3: order cells along the trajectory
HSMM <- orderCells(HSMM)
#HSMM <- orderCells(HSMM, root_state = GM_state(HSMM))
save(HSMM,file=paste0(home,"data_DGE/monocle_2.24GermClusters1kgenes_clustermarkerbymonocle.Robj"))
HSMM2=HSMM

pdf(paste0(file,"2.trajectory_24GermClusters.pdf"),width=6,height=6)
plot_cell_trajectory(HSMM, color_by = "V2")
plot_cell_trajectory(HSMM, color_by = "celltype")
plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()

pData=pData(HSMM)
write.table(pData,paste0(home,"data_DGE/monocle_24GermClusters_pData.txt"),sep="\t",quote=F,row.names=T,col.names=T)

### Finding Genes that Distinguish State
diff_test_res <- differentialGeneTest(HSMM[HSMM_expressed_genes,],
                                             fullModelFormulaStr = '~State',
                                             cores = 8)  
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 5,
                        cores = 8,
                        show_rownames = T)

# used Seurat FindMarkers to search for markers for each of the 5 States
# searched Gene Ontology for 5 states

###### visualize 12 germ cell clusters in Monocle pseudotime and markers - Related to Figure S2C top right panel
### 12 ordered germ cell clusters obtained by Louvain-Jaccard clustering
load(file=paste0(home,"data_DGE/24GermClusters1kgenes_ReScaled_HVG.Robj"))
dge24HVG=dge
which(rownames(pData(HSMM))!=rownames(dge24HVG@data.info))
pData(HSMM)$SPG45new9clusters=dge24HVG@data.info$SPG45new9clusters
file="figAug2017_MouseAdultST24_ReCluster24GermClustersLargeCells1kgenes/Monocle_"
table(round(pData(HSMM)$Pseudotime,0),pData(HSMM)$SPG45new9clusters)
library(RColorBrewer)
myBrewerPalette=brewer.pal(12,"Paired")
# use mannual color by +scale_color_manual(values=myBrewerPalette) in ggplot
pdf(paste0(file,"2.trajectory_Pseudotime.pdf"),height=4.5,width=4)
plot_cell_trajectory(HSMM, color_by = "Pseudotime",show_tree=FALSE,show_branch_points=FALSE,show_backbone=TRUE)
dev.off()
# save as Figure S2C top right panel 
pdf(paste0(file,"2.trajectory_12Clusters_facetwrap.pdf"),width=6)
plot_cell_trajectory(HSMM, color_by = "SPG45new9clusters",show_branch_points=FALSE) + facet_wrap(~SPG45new9clusters, nrow = 3)
dev.off()

### Visualize known markers for germ cells in Monocle
markers=c("Zbtb16","Sall4","Sohlh1","Id4","Gfra1","Uchl1","Kit","Stra8","Prdm9","H2afx","Spo11","Hormad1","Piwil1","Sycp3","Spag6","Mns1","Tbpl1","Acrv1","Tssk1","Tnp1","Prm1","Pgk2","Hspa1l")
HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
        num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                  gene_short_name %in% markers))
cds_subset <- HSMM_filtered[my_genes,]
pdf(paste0(file,"2.12Clusters_knownmarkers.pdf"),width=12)
plot_genes_in_pseudotime(cds_subset, color_by = "SPG45new9clusters",ncol=5)
dev.off()

###### compare pseudotemporal ordering from Waterfall and Monocle with the 12 ordered germ cell (GC) clusters and SOM 20 clusters - Related to Figure S2C bottom panels
setwd("C:/Users/qzm/Desktop/AdultMouseST/plot")
monocle=read.table("data_DGE/Monocle_pseudotime.txt",header=T)
waterfall=read.table("data_DGE/Waterfall_24GermClusters_pseudotime_12clusters.txt")
som=read.table("data_DGE/SOM-20-1_NonSpgGermCells1kgenes_HVG.txt")
gc=read.table("data_DGE/MouseAdultST24_12GermClusters_ClusterID.txt")

gc$cb=rownames(gc)
names(som)[1]=names(monocle)[1]="cb"
waterfall$cb=rownames(waterfall)
names(monocle)[2]="Monocle"
waterfall=waterfall[,-c(1:2)]
names(waterfall)[1]="Waterfall"
som=som[,-c(3:4)]
names(som)[2]="SOM-20-1"
names(gc)[1]="GC1-12"

a=Reduce(function(x,y)merge(x,y,all=TRUE,by="cb"),list(gc,som,waterfall,monocle))
plot(a$SOM,a$Waterfall,col=rgb(1,0,0,0.1),pch=16,ylim=c(0.2,1))
plot(a$GC,a$Waterfall,col=rgb(1,0,0,0.1),pch=16)
# save as Figure S2C bottom left panel
plot(a$SOM,a$Monocle,col=rgb(0,0,1,0.1),pch=16)
plot(a$GC,a$Monocle,col=rgb(0,0,1,0.1),pch=16)
# save as Figure S2C bottom right panel
plot(a$Monocle,a$Waterfall,col=rgb(0,0,0,0.1),pch=16)
plot(a$GC,a$SOM,col=rgb(0,0,0,0.1),pch=16,xlim=c(4,12))

