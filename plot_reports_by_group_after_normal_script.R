#This script, given a table where column one has the sample name (with Sample_) and the other has a group ID (1 to X, with max 8 different groups), will make color coded plots.
#Give the path to this table as the second argument. Name of Results directory is first argument. Species name is third argument.

args <- commandArgs(trailingOnly=T)
Raw_dir=getwd()

Results_dir=args[1]

samples_plus_groups<-read.table(args[2])

species<-args[3]

load(paste(Results_dir,"/Rdata",sep=""))

groups<-as.vector(samples_plus_groups$V2)

group_names<-as.vector(samples_plus_groups$V3)

sample_names<-as.vector(samples_plus_groups$V1)

#indices_for_sorting_sample_names<-order(sample_names)

#sample_names<-sample_names[indices_for_sorting_sample_names]

#groups<-groups[indices_for_sorting_sample_names]

#group_names<-group_names[indices_for_sorting_sample_names]

#GC_matrix<-GC_matrix[,indices_for_sorting_sample_names]
#featureCounts<-featureCounts[,indices_for_sorting_sample_names]
#coverage_matrix<-coverage_matrix[,indices_for_sorting_sample_names]
#insert_matrix<-insert_matrix[,indices_for_sorting_sample_names]
#arld<-arld[,indices_for_sorting_sample_names]

samples<-c()

for(i in 1:length(sample_names))
{
samples<-append(samples,paste("Sample_",sample_names[i],sep=""))
}

library(RColorBrewer)
mycol=brewer.pal(8,"Dark2")

library(pheatmap)

#Now, make plots!

color_vs_group_name<-unique(cbind(groups,group_names))

colors_for_legend<-as.numeric(color_vs_group_name[,1])

groups_for_legend<-as.vector(color_vs_group_name[,2])

mycol_for_lineplot<-c()

for(i in 1:length(sample_names))
{
mycol_for_lineplot<-append(mycol_for_lineplot,mycol[groups[i]])
}

#Gene body coverage plotting

png(paste(Results_dir,"/Picard_genebodycoverage.png",sep=""),width=600,height=800,res=200);par(mgp=c(2,0.5,0))

matplot(0:100,coverage_matrix,type="l",lty=1,main="Gene Body Coverage",xlab="Percentile gene body 5'->3'",ylab="Normalized coverage",cex.lab=0.7,cex.axis=0.5,col=mycol_for_lineplot)

legend("bottom",legend=groups_for_legend,col=mycol[colors_for_legend],lwd=1,bty="n",cex=0.25,pt.lwd=1)

dev.off()

pheatmap(coverage_matrix,cluster_rows=FALSE,cluster_cols=TRUE,cellheight=2,fontsize=3,fontsize_row=2,fontsize_col=2,main="Gene Body Coverage",file=paste(Results_dir,"/Picard_genebodycoverage_heatmap.png",sep=""),cellwidth=2)

#Inner distance plotting

png(paste(Results_dir,"/RSeQC_inner_distance.png",sep=""),width=600,height=800,res=200);par(mgp=c(2,0.5,0))

matplot(mids[1:100],insert_matrix[1:100,],lty=1,type="l",main="Inner Distance",xlab="Inner distance (Read 2 start - Read 1 end)",ylab="Read count/millions reads in sample",cex.lab=0.7,cex.axis=0.5,col=mycol_for_lineplot)

legend("topright",legend=groups_for_legend,col=mycol[colors_for_legend],lwd=1,bty="n",cex=0.25,pt.lwd=1)

dev.off()

pheatmap(insert_matrix[21:100,],cluster_rows=FALSE,cluster_cols=TRUE,cellheight=2,fontsize=3,fontsize_row=2,fontsize_col=2,main="Inner Distance",file=paste(Results_dir,"/RSEQC_inner_dist_heatmap.png",sep=""),cellwidth=2)

#GC content plotting

png(paste(Results_dir,"/RSEQC_GCcontent.png",sep=""),width=600,height=800,res=200);par(mgp=c(2,0.5,0))

matplot(GC_percents,GC_matrix,lty=1,type="l",main="GC Content",xlab="GC Content (%)",ylab="Read count/millions reads in sample",cex.lab=0.7,cex.axis=0.5,col=mycol_for_lineplot)

legend("topleft",legend=groups_for_legend,col=mycol[colors_for_legend],lwd=1,bty="n",cex=0.25,pt.lwd=1)

dev.off()

GC_heatmap_title<-"Normalized read count\n per GC content value (~10-90%)"

pheatmap(GC_matrix,cluster_rows=FALSE,cluster_cols=TRUE,fontsize_row=2,fontsize_col=2,fontsize=3,cellheight=2,main=GC_heatmap_title,file=paste(Results_dir,"/RSEQC_GC_content_heatmap.png",sep=""),cellwidth=2)

#featureCounts density plotting

max_density<-c()

for(i in 1:ncol(featureCounts))
{
        density<-density(log10(featureCounts[,i]))
	max_density<-append(max_density,max(density$y))
}

max_density<-max(max_density)

png(paste(Results_dir,"/featureCounts_density.png",sep=""),width=1200,height=1200,res=200);par(mgp=c(2,0.5,0))

plot(density(log10(featureCounts[,1])),ylim=c(0,max_density),xlim=c(0,6),lwd=1,col=mycol_for_lineplot[1],xlab="log10(count)",ylab="Density",main="Densities for read counts per gene",cex.lab=0.7,cex.axis=0.5)

for(i in 2:ncol(featureCounts))
{
	lines(density(log10(featureCounts[,i])),col=mycol_for_lineplot[i],lwd=1)
}

legend("topright",legend=groups_for_legend,col=mycol[colors_for_legend],lwd=1,bty="n",cex=0.35,pt.lwd=1)

dev.off()

#Get MDS/PCA info.

d=dist(t(arld))
h=hclust(d)

h$labels<-sample_names

cm = cmdscale(d, 2)

#Plot dendrogram and MDS plot.

#png(paste(Results_dir,"/Dendrogram_unsupervised_clustering_from_featureCounts.png",sep=""),width=1200,height=800,res=200)
#par(mgp=c(4,0.5,0))
#par(mar=c(6,5,4,2))
#par(cex=0.3)
#plot(h,main="Dendrogram unsupervised clustering",las=2)
#dev.off()

#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#library("gplots")
#distsRL <- dist(t(assay(rld)))
#mat <- as.matrix(distsRL)
#rownames(mat)<-colnames(mat)<-sample_names
#png(paste(Results_dir,"/Dendrogram_as_heatmap.png",sep=""),width=1200,height=800,res=200)
#heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13),cexRow=0.5,cexCol=0.5,offsetRow=-0.25,offsetCol=-0.25)
#dev.off()

PCA_analysis<-prcomp(arld)
components<-PCA_analysis$rotation

png(paste(Results_dir,"/PCA.png",sep=""),width=1200,height=800,res=200);par(mgp=c(2,0.5,0))
plot(components[,1],components[,2],xlab=paste("PCA1 explains ",percentage_variance_explained_by_PCA_dimension_one,"% of the variance",sep=""),ylab="PCA2",cex.lab=0.7,cex.axis=0.5,main="PCA plot",pch=19,cex=0.1,col=mycol_for_lineplot)
text(components[,1],components[,2],sample_names,offset=-0.001,pos=3,cex=0.3,col=mycol_for_lineplot)
dev.off()

png(paste(Results_dir,"/MDS_plot_normalized.png",sep=""),width=1200,height=800,res=200)
plot(cm[,1],cm[,2],pch=19,cex.lab=0.7,cex.axis=0.5,xlab="Dimension1",ylab="Dimension2",main="MDS plot",cex=0.1,col=mycol_for_lineplot)
text(cm,sample_names,offset=-0.001,pos=3,cex=0.3,col=mycol_for_lineplot)
dev.off()

save.image(paste(Results_dir,"/plot_by_group.Rdata",sep=""))

#save.image(paste(Results_dir,"/Rdata",sep=""))

