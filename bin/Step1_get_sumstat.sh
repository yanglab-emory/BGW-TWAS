#!/usr/bin/bash

################################################################
################################################################
# Step 1: obtain summary statistics (aka Score Statistics)
# Runs single-variant GWAS on the training sample genotypes and available expression data.
# first extracts gene location information for target gene from geneFile.
################################################################
################################################################

# Variable needed for obtaining summary statistics
###
# --BGW_dir : Specify the directory of BGW-TWAS tool
# --wkdir : Specify a working directory
# --gene_name : Specify the gene name that should be the same used in `GeneExpFile`
# --GeneExpFile : Specify gene expression file directory
# --geno_dir : Specify the directory of all genotype files
# --LDdir : Specify the directory of all LD files
# --Genome_Seg_Filehead : Specify the file containing the fileheads of all genome segmentations
# --GTfield : Specify the genotype format in the vcf file that should be used: "GT" (default) or e.g., "DS" for dosage
# --num_cores : Specify the number of parallele sessions

#################################
VARS=`getopt -o "" -a -l \
BGW_dir:,wkdir:,gene_name:,GeneExpFile:,geno_dir:,LDdir:,Genome_Seg_Filehead:,GTfield:,num_cores:,clean_output: \
-- "$@"`


if [ $? != 0 ]
then
    echo "Please provide required input arguments. Terminating....." >&2
    exit 1
fi

eval set -- "$VARS"

while true
do
    case "$1" in
        --BGW_dir|-BGW_dir) BGW_dir=$2; shift 2;;
        --wkdir|-wkdir) wkdir=$2; shift 2;;
        --gene_name|-gene_name) gene_name=$2; shift 2;;
        --GeneExpFile|-GeneExpFile) GeneExpFile=$2; shift 2;;
        --geno_dir|-geno_dir) geno_dir=$2; shift 2;;
        --LDdir|-LDdir) LDdir=$2; shift 2;;
        --Genome_Seg_Filehead|-Genome_Seg_Filehead) Genome_Seg_Filehead=$2; shift 2;;
        --GTfield|-GTfield) GTfield=$2; shift 2;;
        --num_cores|-num_cores) num_cores=$2; shift 2;;
        --clean_output|-clean_output) clean_output=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done

##########################################
# Setting Default Input Argument Values
##########################################
GTfield=${GTfield:-"GT"}
num_cores=${num_cores:-1}
clean_output=${clean_output:-1}

if [ -s ${Genome_Seg_Filehead} ]; then
    num_segments=`wc -l ${Genome_Seg_Filehead} | awk '{print $1}'`
    if [ $num_segments -gt 0 ] ; then
        echo ${gene_name} with ${num_segments} genome blocks
    else
        echo ${Genome_Seg_Filehead} has no genome blocks, one per row. Expecting at least one genome block.
        exit
    fi
else
    echo ${Genome_Seg_Filehead} is empty. Please check.
    exit
fi

echo GTfield = $GTfield , number of cores = $num_cores

#### Create work/output directory if not existed
mkdir -p ${wkdir}
mkdir -p ${LDdir}
echo LD directory ${LDdir}

# Set directory for single variant eQTL summary statistics (score statistics)
Score_dir=${wkdir}/${gene_name}_scores
mkdir -p ${Score_dir}
echo Summary score statistic directory ${Score_dir}
cd ${wkdir}
# echo ${wkdir}

# GeneExpFile columns are Gene Name, Chr, Pos, start, end, expr_data

# the following creates a phenotype file for target gene expression trait that includes subject IDs in the first column and gene expression levels in the second column.
if [ -s ${GeneExpFile} ] ; then
    head -1 ${GeneExpFile} | awk '{$1=$2=$3=$4=$5=""; print substr($0,6)}' | tr ' ' '\n' > temp_ID.txt

    awk -F'[\t]' -v gene=${gene_name} '$5==gene{$1=$2=$3=$4=$5=""; print substr($0,6); exit; }' ${GeneExpFile} | tr ' ' '\n' > exp_temp.txt
    paste temp_ID.txt exp_temp.txt > ${wkdir}/${gene_name}_exp_trait.txt

    pheno=${wkdir}/${gene_name}_exp_trait.txt
else
    echo ${GeneExpFile} is empty. Please provide a valid gene expression file.
fi

# calculate variance of the gene expression trait
pv=$(awk '{delta=$2; sum+=$2; ssq+=(delta - avg/NR)^2} END {print ssq/(NR-1)}' ${wkdir}/${gene_name}_exp_trait.txt)
echo quantitative gene expression trait variance = $pv
echo -e ${gene_name} '\t' ${pv} > ${wkdir}/${gene_name}_geneExp_var.txt
rm -f temp_ID.txt exp_temp.txt

## Run in parallele with specified number of processes by -P
seq 1 ${num_segments}  | xargs -I % -n 1 -P ${num_cores} sh ${BGW_dir}/bin/get_score_stat.sh ${pheno} ${geno_dir} ${Score_dir} ${BGW_dir} ${LDdir} ${Genome_Seg_Filehead} % ${GTfield}

if [ $clean_output -eq 1 ] ; then
rm -fr ${wkdir}/${gene_name}_scores/output
fi

echo Step 1 complete for generating eQTL summary statistics!

exit
