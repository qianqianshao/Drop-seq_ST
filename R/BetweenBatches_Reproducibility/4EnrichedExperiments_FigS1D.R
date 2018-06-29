### R script to compare each enriched/depleted experiment against 6 ST datasets by Qianyi on 12/20/2017, related to Figure S1D
### 1 unbiased Experiment group:
# 0) ST                - 6 datasets
### 4 Depleted/Enriched Experiment groups:
# 1) 1n-Depleted       - 2 datasets
# 2) SpG-enriched      - 3 datasets
# 3) INT-enriched      - 6 datasets
# 4) Sertoli-enriched  - 8 datasets

### load merged 25 ST data
load(file =paste0(home,"data_DGE/MouseAdultST25Somatic.Robj"))
dge25=dge

######## 2D Density plot for each enriched experiment group using PC1-PC2
library(gplots)
### 5 Types of Experiments: 
ST=levels(dge@data.info$orig.ident)[1:6]
depleted1n=levels(dge@data.info$orig.ident)[7:8]
spg=levels(dge@data.info$orig.ident)[9:11]
int=levels(dge@data.info$orig.ident)[12:17]
ser=levels(dge@data.info$orig.ident)[18:25]

sets=list(ST,depleted1n,spg,int,ser)
setslabel=c("ST","1nDepleted","SPG","INT","SER")

numPCs=c(1,2)
range(dge@pca.rot[,numPCs])
range(dge@pca.rot[,numPCs[1]])
range(dge@pca.rot[,numPCs[2]])
xyrange=c(-57,36,-48,67)

numPCs=c(1,3)
range(dge@pca.rot[,numPCs])
range(dge@pca.rot[,numPCs[1]])
range(dge@pca.rot[,numPCs[2]])
xyrange=c(-57,36,-25,59)

numPCs=c(1,2)
xyrange=c(-57,36,-48,67)

##### ST6
h6PCs=dge@pca.rot[grepl(paste(ST,collapse="|"),rownames(dge@pca.rot)),numPCs]
dgefile="figDec2017_MouseAdultST25/All_"

jpeg(file=paste0(dgefile,"hist2d_6ST_samerange_PC",paste(numPCs,collapse="PC"),"raw.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
h6=hist2d_range(h6PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main="6ST",cex.main=2)
dev.off()
# log
h6log=log(h6$counts+1)
jpeg(file=paste0(dgefile,"hist2d_6ST_samerange_PC",paste(numPCs,collapse="PC"),"log.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h6log,col=redblue100,main="",axes=F,an=F,cex.main=2)
dev.off()

write.table(h6$counts,paste0(dgefile,"hist2d_50by50_6ST_samerange_PC",paste(numPCs,collapse="PC"),"raw.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(h6log,paste0(dgefile,"hist2d_50by50_6ST_samerange_PC",paste(numPCs,collapse="PC"),"log.txt"),row.names=T,col.names=T,quote=F,sep="\t")


###### heatmap of all datasets in each of the 4 depleted/enriched experiment groups compared against 6 ST batches - Figure S1D
for(i in 2:5){
set=sets[[i]]
setlabel=setslabel[i]
dgefile=paste0("figDec2017_MouseAdultST25/",setlabel,"_")

h2PCs=dge@pca.rot[grepl(paste(set,collapse="|"),rownames(dge@pca.rot)),numPCs]

### plot together
jpeg(file=paste0(dgefile,"hist2d_samerange_PC",paste(numPCs,collapse="PC"),"all.jpeg"),res=300,height=1900,width=650)
par(mfrow=c(3,1),mar=c(0,1,1,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=setlabel,cex.main=2)
# log
h2log=log(h2$counts+1)
image(h2log,col=redblue100,main=setlabel,axes=F,an=F,cex.main=2)
# ratio against 6 ST batches with symmetric color scale
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=setlabel,axes=F,an=F,cex.main=2)
dev.off()

### plot separately
jpeg(file=paste0(dgefile,"hist2d_samerange_PC",paste(numPCs,collapse="PC"),"raw.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=setlabel,cex.main=2)
dev.off()
# log
h2log=log(h2$counts+1)
jpeg(file=paste0(dgefile,"hist2d_samerange_PC",paste(numPCs,collapse="PC"),"log.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2log,col=redblue100,main=setlabel,axes=F,cex.main=2)
dev.off()
# ratio
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
jpeg(file=paste0(dgefile,"hist2d_samerange_PC",paste(numPCs,collapse="PC"),"by6STratio.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=setlabel,axes=F,cex.main=2)
dev.off()

### write table
write.table(h2$counts,paste0(dgefile,"hist2d_50by50_",setlabel,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(h2by6STratio,paste0(dgefile,"hist2d_50by50_",setlabel,"_samerange_PC",paste(numPCs,collapse="PC"),"logratio.txt"),row.names=T,col.names=T,quote=F,sep="\t")
}
### save as Figure S1D

###### evaluate each individual batch from each experiment group
h2loglist=h2by6STratiolist=zlist=list()
for(j in 1:5){
set=sets[[j]]
setlabel=setslabel[j]
dgefile=paste0("figDec2017_MouseAdultST25/",setlabel,"_")

for(i in 1:length(set)){
seti=set[i]
h2PCs=dge@pca.rot[which(gsub("_.*","",rownames(dge@pca.rot))==seti),numPCs]
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=seti,cex.main=2)
# calculate log
h2log=log(h2$counts+1)
# calculate ratio against 6 ST
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
h2loglist[[i]]=h2log
h2by6STratiolist[[i]]=h2by6STratio
zlist[[i]]=z
}

if(j==1){
#plot log-transformed heatmap for each of 6 ST batches
jpeg(file=paste0(dgefile,"hist2d_samerange_PC",paste(numPCs,collapse="PC"),"all.jpeg"),res=300,height=1200*length(set),width=1100)
par(mfrow=c(length(set),1),mar=c(0,1,2,1),mgp=c(1.5,0.5,0))
# log
for(i in 1:length(set)){
	h2log=h2loglist[[i]]
	z=zlist[[i]]
	image(h2log,col=redblue100,main=seti,axes=F,an=F,cex.main=2)
}
dev.off()
} else {
#plot ratio against 6ST for each depleted/enriched batch
jpeg(file=paste0(dgefile,"hist2d_samerange_PC",paste(numPCs,collapse="PC"),"all.jpeg"),res=300,height=1200*length(set),width=1100)
par(mfrow=c(length(set),1),mar=c(0,1,2,1),mgp=c(1.5,0.5,0))
# ratio
for(i in 1:length(set)){
	h2by6STratio=h2by6STratiolist[[i]]
	z=zlist[[i]]
	image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=seti,axes=F,an=F,cex.main=2)
}
dev.off()
}
}
