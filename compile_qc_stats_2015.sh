#!/bin/bash

############################
#
#For RNA-Seq samples, make "QC_statistics.txt" file with relevant information.
#Mapping rates + assorted other QC info
# New York Genome Center
# Heather Geiger (hmgeiger@nygenome.org)
#Version 04-09-2015
############################

#Accept as arguments the name of the results directory, a list of samples (text file, sample names without Sample_), and whether standard inside East River, standard outside East River,  or PDX.

if [ $# -eq 0 ]
then
echo "Here are the arguments for this script."
echo "All options are required."
echo "-s sample_list.txt (Text file with list of samples, one sample per line, no Sample_ in front of names)"
echo "-r Results_directory_name (Directory will be created under the current working directory)"
echo "-d Data source. Acceptable options are ER (from EastRiver), external (outside East River, not PDX), or PDX."
exit 0
fi

#Parse options. Output header for file.

while getopts s:r:d: option
do
case "${option}"
in
s) sample_list=${OPTARG};;
r) Results_dir=${OPTARG};;
d) data_source=${OPTARG};;
esac
done

seq_type=$data_source

Raw_dir=`pwd`

if [ ! -e $Results_dir ];then mkdir $Results_dir;fi

if [ $data_source != PDX ]
then
echo "Tot_reads(M) rRNA(%) Map(%) UQ_map(%) Gene_assn(%) Strand(%) 5prime_cov 3prime_cov Mean_GC(%) Dup(%) Mean_inner_dist CDS(%) UTR(%) Intronic(%) Intergenic(%)" | awk '{ OFS="\t"}{$1=$1;print $0}' > ${Results_dir}/QC_statistics.txt
else
echo "Tot_reads(M) rRNA(%) Map_after_clean(%) UQ_map_after_clean(%) Gene_assn(%) Strand(%) 5prime_cov 3prime_cov Mean_GC(%) Dup(%) Mean_inner_dist CDS(%) UTR(%) Intronic(%) Intergenic(%) Map_unique_mouse(%)" | awk '{ OFS="\t"}{$1=$1;print $0}' > ${Results_dir}/QC_statistics.txt
fi

#Since East River outputs by lane, rather than one total file, we are going to combine these files and output one file for later on.

if [ $data_source = ER ]
then

for name in `cat $sample_list`

do

sample=Sample_${name}

#For rRNA logs, if there is more than one file, take the last line from all files and output to a new file.
#East River has also been varying in how it names the files. So look for both options.

rRNA_lines_type_one=`cat ${sample}/align/analysis/bowtie2/*/*.log | wc -l`

rRNA_lines_type_two=`cat ${sample}/align/analysis/bowtie2/*.log | wc -l`

if [ ! -e ${sample}/Stats ];then mkdir ${sample}/Stats;fi

if [ $rRNA_lines_type_one -gt 0 ]
then
cat ${sample}/align/analysis/bowtie2/*/*.log | grep overall > ${sample}/Stats/${sample}_rRNA_mapping.txt
elif [ $rRNA_lines_type_two -gt 0 ]
then
cat ${sample}/align/analysis/bowtie2/*.log | grep overall > ${sample}/Stats/${sample}_rRNA_mapping.txt
fi

#Cat all the STAR logs. Later can parse them based on the descriptions for the different lines.

cat ${sample}/align/analysis/star/*STAR/Log.final.out > ${sample}/Stats/${sample}_STAR_logs_concatenated.txt

done

fi

for name in `cat $sample_list`
do

#Calculate rRNA the same way regardless of data source. Take all the lines with format "X% overall alignment rate", then average them.
#If file does not exist, just skip looking at this statistic for now.

sample=Sample_${name}

if [ -s ${sample}/Stats/${sample}_rRNA_mapping.txt ]
then
rRNA_rate=`grep overall ${sample}/Stats/${sample}_rRNA_mapping.txt | awk '{sum+=$1} END {printf("%.2f",sum/NR)}'`
else
rRNA_rate="NA"
fi

#For non-PDX, calculate mapping rates the same way.
#Parse the STAR log(s) to get input reads and uniquely and multiply mapped reads, then calculate rates.

if [ $seq_type = ER ];then STAR_file=${sample}/Stats/${sample}_STAR_logs_concatenated.txt;fi
if [ $seq_type = external ];then STAR_file=${sample}/STAR_alignment/*Log.final.out;fi
if [ $seq_type = PDX ];then STAR_file=${sample}/STAR_alignment_to_human_after_remove_mouse/${sample}_vs_human_after_remove_mouse_reads_Log.final.out;fi

if [ $seq_type != PDX ]
then

total_reads=`awk '$3 ~ /input/' $STAR_file | awk '{sum+=$NF} END {print sum}'`
unique_reads=`awk '$1 ~ /Uniquely/ && $4 ~ /number/' $STAR_file | awk '{sum+=$NF} END {print sum}'`
mapped_reads=`awk '($1 ~ /Uniquely/ && $4 ~ /number/) || ($1 ~ /Number/ && $6 ~ /multiple/)' $STAR_file | awk '{sum+=$NF} END {print sum}'`

mapping_rate=`echo $mapped_reads $total_reads | awk '{printf("%.2f",($1 / $2)*100)}'`
unique_rate=`echo $unique_reads $total_reads | awk '{printf("%.2f",($1 / $2) * 100)}'`

total_reads=`echo $total_reads | awk '{printf("%.2f",$1/1000000)}'`

fi

#For PDX, mapping rate and uniquely mapping rate are after removing mouse reads, so require special calculations. Also need to calculate mouse unique rate.

if [ $seq_type = PDX ]
then

total_reads=`awk '$3 ~ /input/' $STAR_file | awk '{sum+=$NF} END {print sum}'`
unique_reads=`awk '$1 ~ /Uniquely/ && $4 ~ /number/' $STAR_file | awk '{sum+=$NF} END {print sum}'`
mapped_reads=`awk '($1 ~ /Uniquely/ && $4 ~ /number/) || ($1 ~ /Number/ && $6 ~ /multiple/)' $STAR_file | awk '{sum+=$NF} END {print sum}'`

mouse_uniques=`awk '{print $NF}' ${sample}/STAR_alignment_to_mouse_and_human/*counts.txt`

total_minus_mouse=`echo $total_reads $mouse_uniques | awk '{print $1 - $2}'`

mapping_rate=`echo $mapped_reads $total_minus_mouse | awk '{printf("%.2f",($1/$2)*100)}'`
unique_rate=`echo $unique_reads $total_minus_mouse | awk '{printf("%.2f",($1/$2)*100)}'`

mouse_unique_rate=`echo $mouse_uniques $total_reads | awk '{printf("%.2f",($1/$2)*100)}'`

fi

#Dup rate, 5' and 3' bias statistics, GC content, strandedness, gene assignment rate, read distribution, and inner distance are calculate the same no matter what.
#Works for single-end data as well. Inner distance will just be listed as "NA".

if [ $seq_type != ER ]
then
RNAMetrics_file=${sample}/Stats/${sample}_RNAMetrics.txt
dup_file=${sample}/Stats/${sample}_MarkDuplicates.metrics.txt
GC_file=${sample}/Stats/${sample}.GC.xls
inner_dist_file=${sample}/Stats/${sample}.inner_distance_freq.txt
featureCounts_file=${sample}/featureCounts/${sample}_counts.txt.summary
else
RNAMetrics_file=${sample}/qc/${name}.RNAMetrics.metrics
dup_file=${sample}/qc/${name}.dedup.metrics
GC_file=${sample}/qc/${name}.GC.xls
inner_dist_file=${sample}/qc/${name}.inner_distance_freq.txt
featureCounts_file=${sample}/qc/${name}_counts.txt.summary
fi

#Calculate assorted statistics.

duplication_rate=`awk '{if(NR==8)printf("%.2f",$8*100)}' $dup_file`

strandedness=`awk '{if(NR == 8)printf("%.2f",$16*100)}' $RNAMetrics_file`

mean_GC=`awk '{if(NR>1){GC_values_times_reads+=$1*$2;reads+=$2}} END {printf("%.2f",GC_values_times_reads/reads)}' $GC_file`

fiveprime_mean_coverage=`awk '{if(NR > 11 && $1 >= 11 && $1 <= 30)sum+=$2} END {printf("%.2f",sum/20)}' $RNAMetrics_file`

threeprime_mean_coverage=`awk '{if(NR > 11 && $1 >= 71 && $1 <= 90)sum+=$2} END {printf("%.2f",sum/20)}' $RNAMetrics_file`

mean_inner_dist=`awk '{inner_distance_values_times_reads+=(($1 + 2.5)*$3);reads+=$3} END {printf("%.2f",inner_distance_values_times_reads/reads)}' $inner_dist_file`

gene_assignment_rate=`awk '{if(NR == 2)assigned=$NF;if(NR == 2 || NR == 3 || NR == 5)unique_reads+=$NF} END {printf("%.2f",(assigned/unique_reads)*100)}' $featureCounts_file`

CDS_rate=`awk '{if(NR == 8)printf("%.2f",$10*100)}' $RNAMetrics_file`
UTR_rate=`awk '{if(NR == 8)printf("%.2f",$11*100)}' $RNAMetrics_file`
intronic_rate=`awk '{if(NR == 8)printf("%.2f",$12*100)}' $RNAMetrics_file`
intergenic_rate=`awk '{if(NR == 8)printf("%.2f",$13*100)}' $RNAMetrics_file`

#Output to file. Only difference will be to add mouse unique rate for PDX.

if [ $seq_type != PDX ]
then
echo $name $total_reads $rRNA_rate $mapping_rate $unique_rate $gene_assignment_rate $strandedness $fiveprime_mean_coverage $threeprime_mean_coverage $mean_GC $duplication_rate $mean_inner_dist $CDS_rate $UTR_rate $intronic_rate $intergenic_rate | awk '{ OFS="\t"}{$1=$1;print $0}' >> $Results_dir/QC_statistics.txt
else
echo $name $total_reads $rRNA_rate $mapping_rate $unique_rate $gene_assignment_rate $strandedness $fiveprime_mean_coverage $threeprime_mean_coverage $mean_GC $duplication_rate $mean_inner_dist $CDS_rate $UTR_rate $intronic_rate $intergenic_rate $mouse_unique_rate | awk '{ OFS="\t"}{$1=$1;print $0}' >> $Results_dir/QC_statistics.txt
fi

done



