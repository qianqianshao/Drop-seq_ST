### Comparison between Stevant Gonad 2018 data (6 somatic clusters) and our adult testis data (11 cell types) by Qianyi on 6/4/2018
### Related to Figure S5B

### Stevant et al Cell Reports 2018 paper for Fetal Gonad 
# Fluidigm C1 single-cell RNA-seq
Downloaded:
1.  GSE97519: Raw Gene Expression Matrix for All 16,459 genes and 400 Cells
2.  TableS2: 2802 markers, RPKM centroids for 6 clusters 

### Our 11 major cell types for Adult mouse testis
MouseAdultST25_11celltypes_centroid_UMI20cell15genes.txt
data format: log(mean(counts-per-10k)+1), 24947 genes (>20 UMIs and >15 cells), 11 cell type centroids

## path
home="C:/Users/qzm/Desktop/AdultMouseST/DropseqMS/ResponseToReviewers/Reviewer5Comments/ComparisonWithOur25ST_Stevant_Gonad_2018/"
setwd(home)
#source("C:/Users/qzm/Desktop/malek/Scripts_Clustering/Rcode_multiplot.R")
#redblue100<-rgb(read.table('C:/Users/qzm/Desktop/AdultMouseST/plot/data_DGE/redblue100.txt',sep='\t',row.names=1,header=T))

###### load Stevant 2018 Fetal Gonad data
### load raw data from GSE97519
raw=read.table("GSE97519_XY_NR5A1-eGFP_single-cell_fetal_gonads.txt.gz")
dim(raw) # [1] 16459   400

### load 6 cluster centroids from Table S2 of Stevant's 2018 Cell Reports paper
# data format: mean(RPKM), 2802 markers (rows), 6 cluster centroids (columns)
Fetal=read.table("FetalGonad6ClusterCentroids2802Genes_RPKM.txt",header=T)
anyDuplicated(Fetal[,1]) # 0 (no duplicted marker for each cluster)
rownames(Fetal)=Fetal[,1]
Fetal=Fetal[,-1]
### number of markers for each cluster
table(Fetal$ClusterID)
# C1  C2  C3  C4  C5  C6 
#377 602 625 389 329 480
dim(Fetal) # [1] 2802   10
Fetal[1:2,]
              C1        C2         C3        C4 C5         C6 ClusterID
Afap1l1 56.46611 0.5203821 0.07745985 0.1990415  0 0.09825918        C1
Akr1c14 53.25356 0.6311708 0.02809066 0.0000000  0 0.27971679        C1
        num_cells_expressed pval qval
Afap1l1                  14    0    0
Akr1c14                  11    0    0

###### load our 11 adult testis cell type centroids
### load 11 cluster centroids
# data format: log(mean(counts-per-10k)+1), 24947 genes with >20 UMIs and >15 cells (rows), 11 cell type centroids (columns)
adult=read.table("MouseAdultST25_11celltypes_centroid_UMI20cell15genes.txt",header=T)
dim(adult)   # [1] 24947    11
### change data format to mean(counts-per-1M) for our 11 adult cell types
adult<-expm1(adult)*100
adult[1:2,]
              InnateLymphoid Macrophage Endothelial    Myoid      Leydig
0610005C13Rik        0.00000    0.00000     0.00000  0.00000   0.5760422
0610007P14Rik       50.74989   23.59116    24.34993 15.54316 130.6675736
                Sertoli    Unknown Spermatogonia Spermatocyte RoundSpermatid
0610005C13Rik  3.394725  0.1850365     0.5991046   0.08347371     0.04512236
0610007P14Rik 44.127863 55.2197240    93.2512536  17.04023095     1.72576737
              Elongating
0610005C13Rik  0.1793064
0610007P14Rik  3.6137912

###### Match overlapping genes
overlapped=rownames(Fetal[rownames(Fetal) %in% rownames(adult),])
length(overlapped) # [1] 2711
Fetal<-Fetal[overlapped,]
table(Fetal$ClusterID)
# C1  C2  C3  C4  C5  C6 
#369 577 609 373 318 465
adult<-adult[overlapped,]
all(rownames(adult) == rownames(Fetal)) #TRUE
summary(adult)
 InnateLymphoid        Macrophage         Endothelial            Myoid        
 Min.   :    0.000   Min.   :    0.000   Min.   :    0.000   Min.   :   0.00  
 1st Qu.:    0.000   1st Qu.:    2.085   1st Qu.:    3.289   1st Qu.:   0.00  
 Median :    8.232   Median :   15.922   Median :   19.749   Median :  16.54  
 Mean   :   53.908   Mean   :   81.079   Mean   :   83.133   Mean   :  95.71  
 3rd Qu.:   37.910   3rd Qu.:   57.020   3rd Qu.:   60.832   3rd Qu.:  66.39  
 Max.   :11945.546   Max.   :15752.389   Max.   :10758.956   Max.   :7769.51  
     Leydig             Sertoli            Unknown          Spermatogonia     
 Min.   :    0.000   Min.   :   0.000   Min.   :    0.000   Min.   :   0.000  
 1st Qu.:    2.031   1st Qu.:   2.149   1st Qu.:    3.298   1st Qu.:   1.237  
 Median :   13.199   Median :  14.833   Median :   18.868   Median :   7.295  
 Mean   :  100.691   Mean   :  75.243   Mean   :   98.973   Mean   :  43.069  
 3rd Qu.:   46.049   3rd Qu.:  63.107   3rd Qu.:   69.014   3rd Qu.:  38.080  
 Max.   :10723.963   Max.   :5506.717   Max.   :13017.828   Max.   :3162.954  
  Spermatocyte      RoundSpermatid        Elongating      
 Min.   :   0.000   Min.   :   0.0000   Min.   :   0.000  
 1st Qu.:   0.303   1st Qu.:   0.4683   1st Qu.:   0.466  
 Median :   1.607   Median :   2.2643   Median :   1.973  
 Mean   :  33.415   Mean   :  29.8274   Mean   :  24.437  
 3rd Qu.:  14.252   3rd Qu.:  14.2833   3rd Qu.:   7.240  
 Max.   :4234.084   Max.   :2102.3386   Max.   :3585.702
summary(Fetal[,1:6])
       C1                 C2                 C3                 C4           
 Min.   :   0.000   Min.   :   0.000   Min.   :   0.000   Min.   :   0.0000  
 1st Qu.:   0.000   1st Qu.:   1.945   1st Qu.:   2.322   1st Qu.:   0.7078  
 Median :   0.048   Median :   7.519   Median :   8.780   Median :   5.7475  
 Mean   :  29.452   Mean   :  26.345   Mean   :  28.410   Mean   :  31.4309  
 3rd Qu.:  11.520   3rd Qu.:  20.112   3rd Qu.:  23.320   3rd Qu.:  21.3940  
 Max.   :6527.126   Max.   :4514.734   Max.   :2904.171   Max.   :2892.7538  
       C5                 C6          
 Min.   :   0.000   Min.   :   0.000  
 1st Qu.:   0.041   1st Qu.:   1.358  
 Median :   3.355   Median :   6.542  
 Mean   :  36.724   Mean   :  38.678  
 3rd Qu.:  17.973   3rd Qu.:  23.076  
 Max.   :4103.012   Max.   :4550.528

###### calculate correlation between Stevant Fetal gonad and our adult testis cell type centroids using 2,711 overlapped markers
all=cbind(Fetal[,1:6],adult)
all[1:2,]
#              C1        C2         C3        C4 C5         C6 InnateLymphoid
#Afap1l1 56.46611 0.5203821 0.07745985 0.1990415  0 0.09825918              0
#Akr1c14 53.25356 0.6311708 0.02809066 0.0000000  0 0.27971679              0
#        Macrophage Endothelial Myoid   Leydig   Sertoli   Unknown Spermatogonia
#Afap1l1   54.34019   194.40743     0  0.00000 1.3813719  1.748908      2.595386
#Akr1c14   21.13904    91.41801     0 20.12133 0.4281911 48.436606      0.000000
#        Spermatocyte RoundSpermatid Elongating
#Afap1l1   10.1636579      8.4770181  7.9314580
#Akr1c14    0.0355993      0.1447198  0.1499203

rho=cor(all,method="spearman")

summary(c(rho[7:17,1:6]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2008  0.3142  0.3813  0.3750  0.4293  0.5891
rho[7:17,1:6]
                      C1        C2        C3        C4        C5        C6
InnateLymphoid 0.3983885 0.3894934 0.4108183 0.3490766 0.3761274 0.3860890
Macrophage     0.4377431 0.3812286 0.4831656 0.3802705 0.4204218 0.3814149
Endothelial    0.5578050 0.4179171 0.4723889 0.3505712 0.3937107 0.3813269
Myoid          0.3980319 0.3955420 0.4624032 0.3138771 0.3782922 0.3330490
Leydig         0.3751004 0.3999337 0.4820953 0.4786200 0.5835438 0.5211701
Sertoli        0.3168112 0.3223428 0.3085043 0.4864483 0.4307665 0.5891220
Unknown        0.4247893 0.4350134 0.5861623 0.3979417 0.4705832 0.3968925
Spermatogonia  0.3266064 0.4439891 0.3458253 0.3859342 0.3420178 0.4355580
Spermatocyte   0.2320403 0.3128482 0.2577419 0.2706987 0.2680448 0.3172842
RoundSpermatid 0.2056281 0.2552169 0.2007740 0.2492779 0.2405273 0.2814645
Elongating     0.2463613 0.2733788 0.2256118 0.3150384 0.2927609 0.3719213
# save as Figure S5B
