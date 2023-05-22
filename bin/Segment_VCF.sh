#!/usr/bin/bash

vcf_dir=$1
seg_bed=$2
out_dir=$3

# module load tabix/0.2.6

echo Segment VCF using genome block information BED file: $seg_bed

######### Segment VCFs 
## Write results to temp directory
## You may also copy your input data files to temp directory if needed
cat $seg_bed | while read $line ; do
	chr=$(echo ${line} | awk '{print $1}')
	pos_start=$(echo ${line} | awk '{print $2}')
	pos_end=$(echo ${line} | awk '{print $3}')
	seg_name=$(echo ${line} | awk '{print $4}')

	echo $chr $pos_start $pos_end $seg_name
	echo Segmenting VCF File:  ${vcf_dir}/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_${chr}.recalibrated_variants.vcf.gz

## Tabix genotype
	tabix -h ${vcf_dir}/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_${chr}.recalibrated_variants.vcf.gz ${chr}:${pos_start}-${pos_end} >> ${out_dir}/${seg_name}.vcf
	bgzip -f ${out_dir}/${seg_name}.vcf
	tabix -f ${out_dir}/${seg_name}.vcf.gz

done



