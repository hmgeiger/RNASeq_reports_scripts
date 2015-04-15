#!/bin/bash

#Results_dir=args[1]

#samples_plus_groups<-read.table(args[2])

#species<-args[3]

#load(paste(Results_dir,"/Rdata",sep=""))

#groups<-as.vector(samples_plus_groups$V2)

#group_names<-as.vector(samples_plus_groups$V3)

#sample_names<-as.vector(samples_plus_groups$V1)

#samples<-c()

#for(i in 1:length(sample_names))
#{
#samples<-append(samples,paste("Sample_",sample_names[i],sep=""))
#}

if [ $# -ne 2 ]
then
echo "Usage: $0 Results_dir explanation_file.txt"
echo "This script assumes that the main reports script has already been run, with data saved to an Rdata file"
echo "Text file should be in the following format, with one sample per line:"
echo "Sample_name(without Sample_)<tab>Group_number(1 to up to 8)<tab>Group name"
echo "Use list in same order as used to make original reports"
exit 1
fi

Results_dir=$1
group_file=$2

/data/NYGC/Software/R/R-3.0.0/bin/Rscript /data/NYGC/Resources/RNASeqPipelineResources/scripts/plot_reports_by_group_after_normal_script.R $Results_dir $group_file
