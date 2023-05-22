
###### Segment VCF files according to genome block segmentation information provided by lddetect
bed_file=lddetect_1KG_EUR_hg19.bed

### To segment the genome block with chromosome, start position, end position provided in the 6th row of the BED file
chr=$(awk 'NR==6{print $1}' lddetect_1KG_EUR_hg19.bed)
vcf_file=CHR${chr}.vcf.gz

pos_start=$(awk 'NR==6{print $2}' lddetect_1KG_EUR_hg19.bed)
pos_end=$(awk 'NR==6{print $3}' lddetect_1KG_EUR_hg19.bed)
seg_name=GTEx_WGS_CHR${chr}_${pos_start}_${pos_end}

tabix -h ${vcf_file} ${chr}:${pos_start}-${pos_end} > ${seg_name}.vcf
bgzip -f ${seg_name}.vcf
tabix -f ${seg_name}.vcf.gz

