#!/bin/bash

if [ $# -eq 0 ]
then
echo "Here are the arguments for this script."
echo "-s sample_list.txt (File with one sample per line, no "Sample_", optional)"
echo "-p protocol (totalRNA or mRNA- only needed if flag=TRUE)"
echo "-r Results_dir"
echo "-o organism(human or mouse)"
echo "-f flag (TRUE or FALSE- required- do you want samples to be flagged if stats are bad?)"
echo "-l length- Read length to use for GC plotting. If samples have different read lengths, input the longer one, and script will interpolate for the shorter lengths"
echo "-e EastRiver (TRUE or FALSE- Did data come from East River or through manual run of pipeline?)"
echo "-x PDX- Set equal to TRUE to calculate stats based on PDX data"
exit 1
fi

while getopts p:r:o:f:l:e:s:x: option
do
case "${option}"
in
p) protocol=${OPTARG};;
r) Results_dir=${OPTARG};;
o) organism=${OPTARG};;
f) flag=${OPTARG};;
l) length=${OPTARG};;
e) EastRiver=${OPTARG};;
s) sample_list=${OPTARG};;
x) PDX=${OPTARG};;
esac
done

#Required inputs to R script- sample list, flag, read length, EastRiver, species, results directory.
#Optional- Protocol - only needed if flag equal to true.

if [ $flag = TRUE ]
then
if [ -z $protocol ]
then
echo "Need to specify protocol (totalRNA or mRNA) using -p option if flag=TRUE"
exit 1
fi
fi

if [ -z $Results_dir ]
then
echo "Please specify name of Results directory you want to output to using -r option"
exit 1
fi

if [ -z $organism ]
then
echo "Please specify whether human or mouse using -o option"
exit 1
fi

if [ -z $flag ]
then
echo "Please specify TRUE/FALSE for whether or not you want to flag based on quality using -f option"
exit 1
fi

if [ -z $length ]
then
length=100
echo "Warning: Read length of 100 will be used to run interpolation for GC content plotting"
fi

if [ -z $EastRiver ]
then
echo "Please specify TRUE/FALSE for whether or not data comes from EastRiver"
exit 1
fi

if [ $EastRiver = TRUE ]
then
data_type="ER"
fi

if [ $EastRiver = FALSE ]
then
if [ -z $PDX ];then data_type="external";fi
if [ $PDX = TRUE ];then data_type="PDX";fi
fi

#Make results directory, set up HTML stuff.

if [ ! -e $Results_dir ];then mkdir $Results_dir;fi
if [ ! -e $Results_dir/csslib ];then mkdir $Results_dir/csslib;fi
if [ ! -e $Results_dir/jslib ];then mkdir $Results_dir/jslib;fi

if [ -z $sample_list ]
then
ls | grep ^Sample | awk -F "Sample_" '{print $2}' > $Results_dir/samples.txt
sample_list=$Results_dir/samples.txt
fi

Raw_dir=`pwd`

cp /data/NYGC/Resources/RNASeqPipelineResources/jslib/*.js $Raw_dir/$Results_dir/jslib
cp /data/NYGC/Resources/RNASeqPipelineResources/csslib/*.css $Raw_dir/$Results_dir/csslib

#Get QC stats table.

/data/NYGC/Resources/RNASeqPipelineResources/scripts/compile_qc_stats_2015.sh -s $sample_list -r $Results_dir -d $data_type

#Run Rscript to make plots and output HTML page.

#Results_dir=args[1]
#species=args[2]
#flag=args[3]
#EastRiver=args[4]
#read_length=args[5]
#sample_names=read.table(args[6])
#if(length(args) == 7){protocol=args[7]}

if [ -z $protocol ]
then
/data/NYGC/Software/R/R-3.0.0/bin/Rscript /data/NYGC/Resources/RNASeqPipelineResources/scripts/reports_2015.R $Results_dir $organism $flag $EastRiver $length $sample_list NONE
else
/data/NYGC/Software/R/R-3.0.0/bin/Rscript /data/NYGC/Resources/RNASeqPipelineResources/scripts/reports_2015.R $Results_dir $organism $flag $EastRiver $length $sample_list $protocol
fi
