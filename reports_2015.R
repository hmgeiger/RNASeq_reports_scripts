args <- commandArgs(trailingOnly=T)

Results_dir=args[1]
species=args[2]
flag=args[3]
EastRiver=args[4]
read_length=args[5]
read_length=as.numeric(read_length)
sample_names=read.table(args[6])
protocol=args[7]

Raw_dir=getwd()

#sample_names=read.table("Results_Apr13_test/samples.txt")

sample_names<-as.vector(sample_names$V1)
samples<-c()
for(i in 1:length(sample_names))
{
	samples<-append(samples,paste("Sample_",sample_names[i],sep=""))
}

#Make HTML page.

#cp /data/NYGC/Resources/RNASeqPipelineResources/HTML_page_templates/QC_tables_header.html $Raw_dir/$Results_dir/QC.html

system(paste("cp /data/NYGC/Resources/RNASeqPipelineResources/HTML_page_templates/QC_tables_header.html ",Results_dir,"/QC.html",sep=""))

#Load libraries.

.libPaths("/data/NYGC/Resources/RNASeqPipelineResources/R_packages")

library(htmlTable,lib.loc="/data/NYGC/Resources/RNASeqPipelineResources/R_packages")
library(knitr,lib.loc="/data/NYGC/Resources/RNASeqPipelineResources/R_packages")

library(RColorBrewer,lib.loc="/data/NYGC/Resources/RNASeqPipelineResources/R_packages")
mycol=brewer.pal(8,"Dark2")

library(pheatmap,lib.loc="/data/NYGC/Resources/RNASeqPipelineResources/R_packages")

library(Gmisc,lib.loc="/data/NYGC/Resources/RNASeqPipelineResources/R_packages")
library(DESeq2,lib.loc="/data/NYGC/Resources/RNASeqPipelineResources/R_packages")

library("gplots",lib.loc="/data/NYGC/Resources/RNASeqPipelineResources/R_packages")

sink(paste(Results_dir,"/R_session_info.txt",sep=""))
sessionInfo()
sink()

#Output QC statistics to HTML page.

if(flag == "TRUE" & protocol == "mRNA")
{
rRNA_fail=10
rRNA_flag=5
CDS_and_UTR_flag=70
CDS_and_UTR_fail=50

dup_flag=25
dup_fail=30

mapping_flag=85
mapping_fail=80
}

if(flag == "TRUE" & protocol == "totalRNA")
{
rRNA_fail=20
rRNA_flag=10
CDS_and_UTR_flag=60
CDS_and_UTR_fail=40

dup_flag=25
dup_fail=30

mapping_flag=85
mapping_fail=80
}

QC_stats<-read.table(paste(Results_dir,"/QC_statistics.txt",sep=""),header=TRUE,row.names=1,check.names=FALSE,sep="\t")

for(i in 1:ncol(QC_stats))
{
QC_stats[,i]<-as.vector(QC_stats[,i])
}

if(flag == "TRUE")
{
	rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix<-c()

	for(i in 1:nrow(QC_stats))
	{
		rRNA_status<-c()
		map_status<-c()
		dup_status<-c()
		CDS_and_UTR_status<-c()
		if(QC_stats[i,2] > rRNA_fail){rRNA_status="FAIL"}
		if(QC_stats[i,2] <= rRNA_fail & QC_stats[i,2] > rRNA_flag){rRNA_status="FLAG"}
		if(QC_stats[i,2] <= rRNA_flag){rRNA_status="PASS"}
		if(QC_stats[i,3] < mapping_fail){map_status="FAIL"}
		if(QC_stats[i,3] >= mapping_fail & QC_stats[i,3] < mapping_flag){map_status="FLAG"}
		if(QC_stats[i,3] >= mapping_flag){map_status="PASS"}
		if(QC_stats[i,10] > dup_fail){dup_status="FAIL"}
		if(QC_stats[i,10] <= dup_fail & QC_stats[i,10] > dup_flag){dup_status="FLAG"}
		if(QC_stats[i,10] <= dup_flag){dup_status="PASS"}
		if((QC_stats[i,12] + QC_stats[i,13]) < CDS_and_UTR_fail){CDS_and_UTR_status="FAIL"}
		if((QC_stats[i,12] + QC_stats[i,13]) >= CDS_and_UTR_fail & (QC_stats[i,12] + QC_stats[i,13]) < CDS_and_UTR_flag){CDS_and_UTR_status="FLAG"}
		if((QC_stats[i,12] + QC_stats[i,13]) >= CDS_and_UTR_flag){CDS_and_UTR_status="PASS"}
		rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix<-rbind(rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix,c(rRNA_status,map_status,dup_status,CDS_and_UTR_status))
	}

	colnames(rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix)<-c("rRNA_status","Map_status","Dup_status","CDS_and_UTR_status")

	QC_stats<-cbind(QC_stats,rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix)

	colors_for_table<-c()

	for(i in 1:nrow(QC_stats))
	{
		if(rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,1] == "FAIL" | rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,2] == "FAIL" | rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,3] == "FAIL" | rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,4] == "FAIL")
        		{colors_for_table<-append(colors_for_table,"#F08080")}
		if(rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,1] != "FAIL" & rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,2] != "FAIL" & rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,3] != "FAIL" & rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,4] != "FAIL" & (rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,1] == "FLAG" | rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,2] == "FLAG" | rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,3] == "FLAG" | rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,4] == "FLAG"))
        		{colors_for_table<-append(colors_for_table,"#FFD700")}
		if(rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,1] == "PASS" & rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,2] == "PASS" & rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,3] == "PASS" & rRNA_mapping_dup_and_CDS_plus_UTR_status_matrix[i,4] == "PASS")
        	{colors_for_table<-append(colors_for_table,"white")}
	}

	for(i in 1:ncol(QC_stats))
	{
		QC_stats[,i]<-as.vector(QC_stats[,i])
	}

	sink(paste(Results_dir,"/QC.html",sep=""),,append=TRUE)
	normal_print(htmlTable(QC_stats,col.rgroup = colors_for_table,label="myTable",align=rep("c",times=ncol(QC_stats) + 1)),results='asis')
	sink()

	write.table(QC_stats,file=paste(Results_dir,"/QC_statistics_with_PASS_FAIL_FLAG_statuses.txt",sep=""),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

}

if(flag == "FALSE")
{

        sink(paste(Results_dir,"/QC.html",sep=""),,append=TRUE)
        normal_print(htmlTable(QC_stats,label="myTable",align=rep("c",times=ncol(QC_stats) + 1)),results='asis')
        sink()

}

#After this, run system command to output plots to QC page.

#sed -e '1,13d' /data/NYGC/Resources/RNASeqPipelineResources/HTML_page_templates/QC_plots.html >> $Raw_dir/$Results_dir/QC.html

system(paste("sed -e '1,13d' /data/NYGC/Resources/RNASeqPipelineResources/HTML_page_templates/QC_plots.html >> ",Results_dir,"/QC.html",sep=""))

#Make plots.

#Make objects and plot for read alignment and read distribution.

stats<-read.table(paste(Results_dir,"/QC_statistics.txt",sep=""),header=TRUE)
distribution<-stats[,12:15]
stats<-stats[,3:4]
Uniquely_mapped<-stats[,2]
Multiply_mapped<-stats[,1] - stats[,2]
Unmapped<- 100 - stats[,1]

mapping_rates<-cbind(Uniquely_mapped,Multiply_mapped,Unmapped)
rownames(mapping_rates)<-rownames(stats)

mapping_rates<-t(as.matrix(mapping_rates))

png(paste(Results_dir,"/read_alignment.png",sep=""),width=1200,height=800,res=200)
par(mar=c(12,5,4,0));par(xpd=TRUE)
par(mgp=c(2,0.5,0))
layout(rbind(c(1,1,1,1,1,1,1,2)))
barplot(mapping_rates,col=mycol,ylab="Percent of total reads",las=2,main="Read Alignment Statistics",cex.names=0.5)
par(mfg=c(1,2))
legend("topright",rownames(mapping_rates),pch=22,pt.bg=mycol,bty="n",x.intersp=0.15,cex=0.5,pt.cex=1)
dev.off()

colnames(distribution)<-c("CDS","UTR","Intron","Intergenic")
distribution<-t(as.matrix(distribution))

png(paste(Results_dir,"/read_distribution.png",sep=""),width=1200,height=800,res=200)
par(mar=c(12,5,4,0));par(xpd=TRUE)
par(mgp=c(2,0.5,0))
layout(rbind(c(1,1,1,1,1,1,1,2)))
barplot(distribution,col=mycol,ylab="Percent of total assigned tags",las=2,main="Read Distribution Statistics",cex.names=0.5)
par(mfg=c(1,2))
legend("topright",rownames(distribution),pch=22,pt.bg=mycol,bty="n",x.intersp=0.15,cex=0.5,pt.cex=1)
dev.off()

#Get matrices for gene body coverage, insert size, featureCounts, and GC content.

if(EastRiver == "FALSE")
{
	for(i in 1:length(samples))
	{
        	file=read.table(paste(Raw_dir,"/",samples[i],"/featureCounts/Sample_",sample_names[i],"_counts.txt",sep=""),skip=2,row.names=1)
        	if(i == 1){featureCounts<-file[,c(1,5:6)]}
        	else{featureCounts[,(i+2)]<-file[,6]}
	}

	colnames(featureCounts)<-c("Chromosome","Length",sample_names)

	featureCounts<-featureCounts[order(rownames(featureCounts)),]

	write.table(featureCounts,paste(Results_dir,"/featureCounts_count_matrix.txt",sep=""),quote=FALSE,row.names=TRUE,sep="\t")

	featureCounts_with_chromosome_and_length_info<-featureCounts
	featureCounts<-featureCounts[,3:ncol(featureCounts)]

	save.image(paste(Results_dir,"/Rdata",sep=""))

	coverage_matrix<-c()

	for(i in 1:length(samples))
	{
		file=read.table(paste(Raw_dir,"/",samples[i],"/Stats/Sample_",sample_names[i],"_RNAMetrics.txt",sep=""),header=FALSE,sep="\t",skip=11)
		coverage_matrix<-cbind(coverage_matrix,file[,2])
	}

	rownames(coverage_matrix)<-0:100
	colnames(coverage_matrix)<-sample_names

	save.image(paste(Results_dir,"/Rdata",sep=""))

	insert_matrix<-c()

	for(i in 1:length(samples))
	{
        	file=read.table(paste(Raw_dir,"/",samples[i],"/Stats/Sample_",sample_names[i],".inner_distance_freq.txt",sep=""),header=FALSE)
        	insert_matrix<-cbind(insert_matrix,file[,3]/(sum(file[,3])/1000000))
	}

	starts=seq(from=-250,to=495,by=5)
	ends=seq(from=-245,to=500,by=5)
	mids=starts + 2.5

	intervals<-c()

	for(i in 1:length(starts))
	{
        	intervals[i]<-paste(starts[i],"/",ends[i],sep="")
	}

	rownames(insert_matrix)<-intervals
	colnames(insert_matrix)<-sample_names

	GC_matrix<-c()

	for(i in 1:length(samples))
	{
		GC_content=read.table(paste(Raw_dir,"/",samples[i],"/Stats/Sample_",sample_names[i],".GC.xls",sep=""),header=FALSE,skip=1)
		GC_content<-GC_content[order(GC_content[,1]),]
		GC_interpolated<-approx(GC_content$V1,GC_content$V2,seq(from=100,to=read_length*100,by=100)/read_length,rule=2)
		GC_interpolated_normalized_reads<-GC_interpolated$y/(sum(GC_interpolated$y)/1000000)
		GC_matrix<-cbind(GC_matrix,GC_interpolated_normalized_reads)
	}

	GC_percents<-round(seq(from=100,to=read_length*100,by=100)/read_length,digits=2)

	rownames(GC_matrix)<-round(seq(from=100,to=read_length*100,by=100)/read_length,digits=2)
	colnames(GC_matrix)<-sample_names

	save.image(paste(Results_dir,"/Rdata",sep=""))

}

if(EastRiver == "TRUE")
{

	for(i in 1:length(samples))
	{
        	file=read.table(paste(Raw_dir,"/",samples[i],"/qc/",sample_names[i],"_counts.txt",sep=""),skip=2,row.names=1)
        	if(i == 1){featureCounts<-file[,c(1,5:6)]}
        	else{featureCounts[,(i+2)]<-file[,6]}
	}

	colnames(featureCounts)<-c("Chromosome","Length",sample_names)

	featureCounts<-featureCounts[order(rownames(featureCounts)),]

	write.table(featureCounts,paste(Results_dir,"/featureCounts_count_matrix.txt",sep=""),quote=FALSE,row.names=TRUE,sep="\t")

	featureCounts_with_chromosome_and_length_info<-featureCounts
	featureCounts<-featureCounts[,3:ncol(featureCounts)]

	save.image(paste(Results_dir,"/Rdata",sep=""))

	coverage_matrix<-c()

	for(i in 1:length(samples))
	{
		file=read.table(paste(Raw_dir,"/",samples[i],"/qc/",sample_names[i],".RNAMetrics.metrics",sep=""),header=FALSE,sep="\t",skip=11)
		coverage_matrix<-cbind(coverage_matrix,file[,2])
	}

	rownames(coverage_matrix)<-0:100
	colnames(coverage_matrix)<-sample_names

	save.image(paste(Results_dir,"/Rdata",sep=""))

	insert_matrix<-c()

	for(i in 1:length(samples))
	{
        	file=read.table(paste(Raw_dir,"/",samples[i],"/qc/",sample_names[i],".inner_distance_freq.txt",sep=""),header=FALSE)
       		insert_matrix<-cbind(insert_matrix,file[,3]/(sum(file[,3])/1000000))
	}

	starts=seq(from=-250,to=495,by=5)
	ends=seq(from=-245,to=500,by=5)
	mids=starts + 2.5

	intervals<-c()

	for(i in 1:length(starts))
	{
        	intervals[i]<-paste(starts[i],"/",ends[i],sep="")
	}

	rownames(insert_matrix)<-intervals
	colnames(insert_matrix)<-sample_names

	GC_matrix<-c()

	for(i in 1:length(samples))
	{
		GC_content=read.table(paste(Raw_dir,"/",samples[i],"/qc/",sample_names[i],".GC.xls",sep=""),header=FALSE,skip=1)
		GC_content<-GC_content[order(GC_content[,1]),]
		GC_interpolated<-approx(GC_content$V1,GC_content$V2,seq(from=100,to=read_length*100,by=100)/read_length,rule=2)
		GC_interpolated_normalized_reads<-GC_interpolated$y/(sum(GC_interpolated$y)/1000000)
		GC_matrix<-cbind(GC_matrix,GC_interpolated_normalized_reads)
	}
	
	GC_percents<-round(seq(from=100,to=read_length*100,by=100)/read_length,digits=2)
	rownames(GC_matrix)<-round(seq(from=100,to=read_length*100,by=100)/read_length,digits=2)
	colnames(GC_matrix)<-sample_names

	save.image(paste(Results_dir,"/Rdata",sep=""))
}

#Plot gene body coverage, insert sizes, GC content, and counts density.

png(paste(Results_dir,"/Picard_genebodycoverage.png",sep=""),width=600,height=800,res=200);par(mgp=c(2,0.5,0))

matplot(0:100,coverage_matrix,col=mycol,lty=1,type="l",main="Gene Body Coverage",xlab="Percentile gene body 5'->3'",ylab="Normalized coverage",cex.lab=0.7,cex.axis=0.5)

if(length(sample_names) <= 16)
{
legend("bottom",legend=sample_names,col=mycol,lwd=1,bty="n",cex=0.25,pt.lwd=1,ncol=2)
}
if(length(sample_names) > 16)
{
text(40,0.4,labels="See heatmap-\ntoo many samples to display legend on plot.",cex=0.3)
}
dev.off()

pheatmap(coverage_matrix,cluster_rows=FALSE,cluster_cols=TRUE,cellheight=2,fontsize=3,fontsize_row=2,fontsize_col=2,main="Gene Body Coverage",file=paste(Results_dir,"/Picard_genebodycoverage_heatmap.png",sep=""),cellwidth=2)

png(paste(Results_dir,"/RSeQC_inner_distance.png",sep=""),width=600,height=800,res=200);par(mgp=c(2,0.5,0))

matplot(mids[1:100],insert_matrix[1:100,],col=mycol,lty=1,type="l",main="Inner Distance",xlab="Inner distance (Read 2 start - Read 1 end)",ylab="Read count/millions reads in sample",cex.lab=0.7,cex.axis=0.5)
if(length(sample_names) <= 16)
{
legend("topright",legend=sample_names,col=mycol,lwd=1,bty="n",cex=0.25,pt.lwd=1,ncol=2)
}
dev.off()

pheatmap(insert_matrix[21:100,],cluster_rows=FALSE,cluster_cols=TRUE,cellheight=2,fontsize=3,fontsize_row=2,fontsize_col=2,main="Inner Distance",file=paste(Results_dir,"/RSEQC_inner_dist_heatmap.png",sep=""),cellwidth=2)

png(paste(Results_dir,"/RSEQC_GCcontent.png",sep=""),width=600,height=800,res=200);par(mgp=c(2,0.5,0))
matplot(GC_percents,GC_matrix,lty=1,type="l",main="GC Content",xlab="GC Content (%)",ylab="Read count/millions reads in sample",cex.lab=0.7,cex.axis=0.5,col=mycol)

if(length(sample_names) <= 16)
{
legend("topleft",legend=sample_names,col=mycol,lwd=1,bty="n",cex=0.25,pt.lwd=1,ncol=2)
}
dev.off()

pheatmap(GC_matrix,cluster_rows=FALSE,cluster_cols=TRUE,fontsize_row=2,fontsize_col=2,fontsize=3,cellheight=2,main="Normalized read count\n per GC content value (~10-90%)",file=paste(Results_dir,"/RSEQC_GC_content_heatmap.png",sep=""),cellwidth=2)

max_density<-c()

for(i in 1:ncol(featureCounts))
{
        density<-density(log10(featureCounts[,i]))
        max_density<-append(max_density,max(density$y))
}

max_density<-max(max_density)

png(paste(Results_dir,"/featureCounts_density.png",sep=""),width=1200,height=1200,res=200);par(mgp=c(2,0.5,0))

plot(density(log10(featureCounts[,1])),ylim=c(0,max_density),xlim=c(0,6),lwd=1,col=mycol[1],xlab="log10(count)",ylab="Density",main="Densities for read counts per gene",cex.lab=0.7,cex.axis=0.5)

for(i in 2:ncol(featureCounts))
{
        lines(density(log10(featureCounts[,i])),col=mycol[((i-1)%%8)+1],lwd=1)
}

if(length(sample_names) <= 16)
{
legend("topright",legend=sample_names,col=mycol,lwd=1,bty="n",cex=0.35,pt.lwd=1,ncol=2)
}
dev.off()

#Plot cumulative curve and heatmap.

cumulative_percents_matrix_for_featureCounts<-c()

for(i in 1:ncol(featureCounts))
{
        cumulative_sums<-cumsum(sort(featureCounts[,i],decreasing=T))*100
        cumulative_percents<-cumulative_sums/sum(featureCounts[,i])
        cumulative_percents_matrix_for_featureCounts<-cbind(cumulative_percents_matrix_for_featureCounts,cumulative_percents)
}

save.image(paste(Results_dir,"/Rdata",sep=""))

png(paste(Results_dir,"/Cumulative_curve.png",sep=""),width=600,height=800,res=200);par(mgp=c(2,0.5,0))
matplot(cumulative_percents_matrix_for_featureCounts,log="x",type="l",lty=1,col=mycol,ylim=c(0,100),lwd=0.5,main="Cumulative curve",xlab="rank (log10)",ylab="% of the library",cex.lab=0.7,cex.axis=0.5)
if(length(samples) < 25)
{
legend("topleft",legend=sample_names,col=mycol,lwd=1,bty="n",cex=0.25,pt.lwd=1,ncol=2)
}
dev.off()

colnames(cumulative_percents_matrix_for_featureCounts)<-sample_names
rownames(cumulative_percents_matrix_for_featureCounts)<-1:nrow(featureCounts)

pheatmap(cumulative_percents_matrix_for_featureCounts[1:10,],cluster_rows=FALSE,cluster_cols=TRUE,cellheight=2,fontsize=3,fontsize_row=2,fontsize_col=2,main="Cumulative percents",cellwidth=2,file=paste(Results_dir,"/cumulative_percents_from_featureCounts_heatmap.png",sep=""))

#Do species-specific plotting and tables (chromosome/top genes, Xist vs. chrY, geneInfo).

if(species == "human")
{
geneInfo=read.table("/data/NYGC/Resources/ENCODE/Gencode/gencode.v18.annotation_geneInfo.txt",h=T,sep="\t",as.is=T,quote="")
}
if(species == "mouse")
{
geneInfo=read.table("/data/NYGC/Resources/RNASeqPipelineResources/gencode.vM2.annotation.geneInfo.withDescriptions.txt",h=T,sep="\t",as.is=T,quote="\"")
}

rownames(geneInfo)=geneInfo$Ensembl_Gene_ID

featureCounts<-featureCounts[order(rownames(featureCounts)),]
geneInfo<-geneInfo[order(rownames(geneInfo)),]

write.table(geneInfo,paste(Results_dir,"/gene_annotation_info.txt",sep=""),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")

#Get percent of reads per chromosome.

data_by_chr<-cbind(geneInfo$Chromosome,featureCounts)
colnames(data_by_chr)<-c("Chromosome",sample_names)
write.table(data_by_chr,paste(Results_dir,"/temp_chr_data.txt",sep=""),row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
data_by_chr<-read.table(paste(Results_dir,"/temp_chr_data.txt",sep=""),header=TRUE,row.names=1)
data_by_chr<-aggregate(data_by_chr[,-1],data_by_chr["Chromosome"],sum)
data_by_chr2<-data_by_chr[,-1]
rownames(data_by_chr2)<-data_by_chr[,1]
data_by_chr<-data_by_chr2
rm(data_by_chr2)

if(species == "human")
{
data_by_chr<-data_by_chr[c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrM","chrX","chrY"),]
}

if(species == "mouse")
{
data_by_chr<-data_by_chr[c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrM","chrX","chrY"),]
}

data_by_chr_percents<-prop.table(as.matrix(data_by_chr),margin=2)*100

numeric_chr<-colSums(data_by_chr_percents[1:(nrow(data_by_chr) - 3),])
other_chr<-data_by_chr_percents[(nrow(data_by_chr) - 2):nrow(data_by_chr),]
data_by_chr_percents_consolidated<-rbind(numeric_chr,other_chr)

rownames(data_by_chr_percents_consolidated)<-c("chr1-22",rownames(other_chr))
colnames(data_by_chr_percents_consolidated)<-sample_names

write.table(t(data_by_chr_percents),file=paste(Results_dir,"/percent_reads_per_chr.txt",sep=""),quote=FALSE,col.names=TRUE,row.names=sample_names,sep="\t")
write.table(t(data_by_chr_percents_consolidated),file=paste(Results_dir,"/percent_reads_per_chr_consolidated.txt",sep=""),quote=FALSE,col.names=TRUE,row.names=sample_names,sep="\t")

save.image(paste(Results_dir,"/Rdata",sep=""))

#Plot Xist vs. chrY.

Xist<-featureCounts[grep('Xist$',geneInfo$Gene_name,ignore.case=TRUE),]

chrY<-grep('chrY',geneInfo$Chromosome)
coding_genes<-grep('protein_coding',geneInfo$Gene_type)
chrY<-intersect(chrY,coding_genes)
chrY<-colSums(featureCounts[chrY,])

Xist<-as.numeric(Xist)
chrY<-as.numeric(chrY)

png(paste(Results_dir,"/Xist_vs_chrY.png",sep=""),width=1200,height=1200,res=200);par(mgp=c(2,0.5,0))
plot(Xist,chrY,type="p",main="Xist vs. chrY expression",xlab="Xist",ylab="chrY",cex.lab=0.7,cex.axis=0.5,pch=19,cex=0.1,xlim=c(0,max(c(Xist,chrY))),ylim=c(0,max(c(Xist,chrY))))
text(Xist,chrY,sample_names,offset=-0.001,pos=3,cex=0.3)
dev.off()

#Get top genes table and print to results folder.

count<-featureCounts

count<-count[order(rownames(count)),]

m=unique(c(apply(count,2,function(x) order(x,decreasing=T)[1:5])))

library_proportions_per_top_gene<-cbind(geneInfo$Gene_name[m],signif(t(t(count[m,]*100)/apply(count,2,sum)),2))

top_genes_geneIDs<-rownames(library_proportions_per_top_gene)

rownames(library_proportions_per_top_gene)<-library_proportions_per_top_gene[,1]

library_proportions_per_top_gene<-library_proportions_per_top_gene[,2:ncol(library_proportions_per_top_gene)]

save.image(paste(Results_dir,"/Rdata",sep=""))

write.table(t(library_proportions_per_top_gene),paste(Results_dir,"/library_proportions_per_gene_in_top_five.txt",sep=""),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")

geneInfo_top_genes<-geneInfo[top_genes_geneIDs,]
rownames(geneInfo_top_genes)<-geneInfo_top_genes$Gene_name

#Get hemoglobin percentages per sample.

hemoglobin_genes<-grep('globin',geneInfo$description)
hemoglobin_genes<-rownames(geneInfo[hemoglobin_genes,])

hemoglobin_counts<-colSums(featureCounts[hemoglobin_genes,])
hemoglobin_percentages<-round((hemoglobin_counts/colSums(featureCounts))*100,digits=2)


colnames(data_by_chr_percents)<-sample_names
colnames(data_by_chr_percents_consolidated)<-sample_names

data_by_chr_percents<-round(data_by_chr_percents,digits=2)
data_by_chr_percents_consolidated<-round(data_by_chr_percents_consolidated,digits=2)

sink(paste(Results_dir,"/QC.html",sep=""),append=TRUE)
write("<br>",stdout())
write("Proportion of library assigned to genes, by chromosome",stdout())
write("<br>",stdout())
htmlTable(t(as.matrix(data_by_chr_percents)),title="",label="myTable2",align=rep("c",times=ncol(data_by_chr_percents) + 1))
write("<br>",stdout())
htmlTable(t(as.matrix(data_by_chr_percents_consolidated)),title="",label="myTable3",align=rep("c",times=ncol(data_by_chr_percents_consolidated) + 1))
sink()

library_proportions_per_top_gene<-read.table(paste(Results_dir,"/library_proportions_per_gene_in_top_five.txt",sep=""),header=TRUE,row.names=1)
library_proportions_per_top_gene<-t(library_proportions_per_top_gene)
library_proportions_per_top_gene<-rbind(hemoglobin_percentages,library_proportions_per_top_gene)
rownames(library_proportions_per_top_gene)<-c("Globin(series of genes)",rownames(library_proportions_per_top_gene)[2:length(rownames(library_proportions_per_top_gene))])

sink(paste(Results_dir,"/QC.html",sep=""),append=TRUE)
write("<br>",stdout())
write("Proportion of library assigned to genes, per sample per gene, for genes in top 5 for one or more samples, plus globin",stdout())
write("<br>",stdout())
htmlTable(as.matrix(library_proportions_per_top_gene),title="",label="myTable4",align=rep("c",times=nrow(library_proportions_per_top_gene)+1))
sink()

sink(paste(Results_dir,"/QC.html",sep=""),append=TRUE)
write("<br>",stdout())
htmlTable(t(as.matrix(library_proportions_per_top_gene)),title="",label="myTable5",align=rep("c",times=ncol(library_proportions_per_top_gene)+1))
sink()

sink(paste(Results_dir,"/QC.html",sep=""),append=TRUE)
write("<br>",stdout())
write("More info on top expressed genes",stdout())
write("<br>",stdout())
htmlTable(geneInfo_top_genes,title="",label="myTable6",align=rep("c",times=ncol(geneInfo_top_genes)+1))
sink()

#Run system command to fix QC table borders.

#/data/NYGC/Resources/RNASeqPipelineResources/scripts/fix_QC_table_borders.pl $Raw_dir/$Results_dir/QC.html > $Raw_dir/$Results_dir/QC_tables.html_new

#mv $Raw_dir/$Results_dir/QC_tables.html_new $Raw_dir/$Results_dir/QC.html

system(paste("/data/NYGC/Resources/RNASeqPipelineResources/scripts/fix_QC_table_borders.pl ",Results_dir,"/QC.html > ",Results_dir,"/QC_tables.html_new",sep=""))

system(paste("mv ",Results_dir,"/QC_tables.html_new ",Results_dir,"/QC.html",sep=""))
#Normalize counts using DESeq2. Then use normalized counts to make dendrogram/MDS/PCA.
#Note: This section can be commented out to shorten runtime to a minimal level.

count<-featureCounts

design=data.frame(samples=factor(names(count)))
cds <- DESeqDataSetFromMatrix(countData=count,colData=design,formula(~ samples))
cds<-estimateSizeFactors(cds)
normalized_featureCounts<-counts(cds,normalized=T)

rld <- rlogTransformation( cds,blind=TRUE )
arld=assay(rld) 

save.image(paste(Results_dir,"/Rdata",sep=""))

normalized_featureCounts<-normalized_featureCounts[order(rownames(normalized_featureCounts)),]
arld<-arld[order(rownames(arld)),]

colnames(normalized_featureCounts)<-sample_names
colnames(arld)<-sample_names

write.table(normalized_featureCounts,paste(Results_dir,"/DESeq2_normalized_count_matrix.txt",sep=""),quote=FALSE,row.names=TRUE,sep="\t",col.names=sample_names)
write.table(arld,paste(Results_dir,"/DESeq2_regularized_log_transformed.txt",sep=""),quote=FALSE,row.names=TRUE,sep="\t",col.names=sample_names)

#Plot dendrogram and MDS plot.

#Get dendrogram and MDS and PCA info.

d=dist(t(arld))
h=hclust(d)

h$labels<-sample_names

cm = cmdscale(d, 2)

components<-prcomp(arld)
PCA_analysis<-components
components<-components$rotation

percentage_variance_explained_by_PCA_dimension_one<-round(summary(PCA_analysis)$importance[2,1]*100,digits=2)
percentage_variance_explained_by_PCA_dimension_two<-round(summary(PCA_analysis)$importance[2,2]*100,digits=2)

save.image(paste(Results_dir,"/Rdata",sep=""))

#Plot.

png(paste(Results_dir,"/Dendrogram_unsupervised_clustering_from_featureCounts.png",sep=""),width=1200,height=800,res=200)
par(mgp=c(4,0.5,0))
par(mar=c(6,5,4,2))
par(cex=0.3)
plot(h,main="Dendrogram unsupervised clustering",las=2)
dev.off()

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat)<-colnames(mat)<-sample_names
png(paste(Results_dir,"/Dendrogram_as_heatmap.png",sep=""),width=1200,height=800,res=200)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13),cexRow=0.5,cexCol=0.5,offsetRow=-0.25,offsetCol=-0.25)
dev.off()

png(paste(Results_dir,"/PCA.png",sep=""),width=1200,height=800,res=200);par(mgp=c(2,0.5,0))
plot(components[,1],components[,2],xlab=paste("PCA1, explains ",percentage_variance_explained_by_PCA_dimension_one,"% of the variance",sep=""),ylab=paste("PCA2, explains ",percentage_variance_explained_by_PCA_dimension_two,"% of the variance",sep=""),cex.lab=0.7,cex.axis=0.5,main="PCA plot",pch=19,cex=0.1)
text(components[,1],components[,2],sample_names,offset=-0.001,pos=3,cex=0.3)
dev.off()

png(paste(Results_dir,"/MDS_plot_normalized.png",sep=""),width=1200,height=800,res=200)
plot(cm[,1],cm[,2],pch=19,cex.lab=0.7,cex.axis=0.5,xlab="Dimension1",ylab="Dimension2",main="MDS plot",cex=0.1)
text(cm,sample_names,pos=3,offset=-0.001,cex=0.3)
dev.off()

