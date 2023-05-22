#!/usr/bin/bash

# Tool Tabix is required for running this script
# R library data.table and dplyr are required

#################################
VARS=`getopt -o "" -a -l \
BGW_dir:,wkdir:,gene_name:,BGW_weight:,test_geno_dir:,test_geno_filehead:,GTfield:,test_pheno:,num_cores:,clean_output: \
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
		--BGW_weight|-BGW_weight) BGW_weight=$2; shift 2;;
        --test_geno_dir|-test_geno_dir) test_geno_dir=$2; shift 2;;
		--test_geno_filehead|-test_geno_filehead) test_geno_filehead=$2; shift 2;;
		--GTfield|-GTfield) GTfield=$2; shift 2;;
		--test_pheno|-test_pheno) test_pheno=$2; shift 2;;
        --num_cores|-num_cores) num_cores=$2; shift 2;;
        --clean_output|-clean_output) clean_output=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done

##########################################
# Setting Default Input Argument Values
##########################################
GTfield=${Nmcmc:-"GT"}
num_cores=${num_cores:-1}
clean_output=${clean_output:-1}

mkdir -p ${wkdir}/${gene_name}_GReX
cd ${wkdir}/${gene_name}_GReX

## Check if $BGW_weight is generated and there is valid eQTL
if [ ! -s ${BGW_weight} ] ; then
    echo ${BGW_weight} is not generated. Please check Step 3.
    exit
else
    n_eqtl=`grep -v "#" ${BGW_weight} | wc -l  | awk '{print $1}'`
    if [ $n_eqtl -gt 0 ] ; then
        echo Extract genotype vcf file for ${n_eqtl} SNPs in ${BGW_weight} ...
    else
        echo There is no SNPs with non-zero eQTL effect sizes in ${BGW_weight}.
        exit
    fi
fi

## make header row for vcf file from one of the existing genotype files
file_temp=$(head -n1 ${test_geno_filehead})
tabix ${test_geno_dir}/${file_temp}.vcf.gz -H | tail -n1 > ${wkdir}/${gene_name}_GReX/${gene_name}_grex.vcf

# Loop through SNP ID and use tabix
tail -n+2 ${BGW_weight} | while read line ; do
chr=$(echo $line | awk -F'[\t ]' '{print $1}' )
pos=$(echo $line | awk -F'[\t ]' '{print $2}' )
ref=$(echo $line | awk -F'[\t ]' '{print $4}' )
alt=$(echo $line | awk -F'[\t ]' '{print $5}' )
#echo ${chr}:${pos}-${pos} $ref $alt

n_file=$(grep _CHR${chr}_ $test_geno_filehead | wc -l | awk '{print $1}')

if [ $n_file -eq 1 ] ; then
#echo ${file_temp}.vcf.gz
file_temp=$(grep _CHR${chr}_ $test_geno_filehead)
tabix ${test_geno_dir}/${file_temp}.vcf.gz ${chr}:${pos}-${pos} | awk -v ref=${ref} -v alt=${alt} -F"\t" '($4==ref && $5==alt) || ($4==alt && $5==ref) {print }' >> ${wkdir}/${gene_name}_GReX/${gene_name}_grex.vcf
else
    echo There is no file head or multiple file heads with the pattern of _CHR${chr}_ in $test_geno_filehead.
    echo $file_temp
    echo Please combine all genotype in the same chromorome into one VCF file with the pattern of \"_CHR${chr}_\" in the filename and update it in $test_geno_filehead.
fi
done

echo ${wkdir}/${gene_name}_GReX/${gene_name}_grex.vcf is created.

## Generate genotype file (dosages) with ${gene_name}_grex.vcf
nsnp=`grep -v "#" ${wkdir}/${gene_name}_GReX/${gene_name}_grex.vcf | wc -l  | awk '{print $1}'`
if [ $nsnp -gt 0 ] ; then
    ${BGW_dir}/bin/Estep_mcmc -vcf ${wkdir}/${gene_name}_GReX/${gene_name}_grex.vcf -p ${test_pheno} -o ${gene_name}_grex -GTfield ${GTfield} -saveGeno -maf 0
    if [ -s ./output/${gene_name}_grex.geno ] ; then
        # echo Converting test genotype VCF file ${wkdir}/${gene_name}_GReX/${gene_name}_grex.vcf to genotype dosage file ${wkdir}/${gene_name}_GReX/output/${gene_name}_grex.geno is success.
        sed -i 's/#//g' ./output/${gene_name}_grex.geno
        mv -f ./output/${gene_name}_grex.geno ${wkdir}/
        echo Converting test genotype VCF file to genotype dosage file ${wkdir}/${gene_name}_grex.geno is success.
    else
        echo ${wkdir}/${gene_name}_GReX/output/${gene_name}_grex.geno failed to be generated.
        exit
    fi
else
    echo There is no test SNPs with non-zero eQTL weights for gene $gene_name. Please check if you BGW weight file and test VCF files.
    exit
fi

## Remove the "#" from the header row
sed -i 's/#//' $BGW_weight
if [ -s ${wkdir}/${gene_name}_grex.geno ] ; then
    echo Calculating predicted GReX ...
    Rscript ${BGW_dir}/bin/compute_grex.R ${gene_name} ${wkdir}
    echo Test GReX file is generated under ${wkdir}
else
	echo Test genotype dosage file ${wkdir}/${gene_name}_grex.geno file failed to be generated. Please check your input files for ${BGW_dir}/bin/Step4_get_test_grex.sh.
	exit
fi

if [ $clean_output -eq 1 ] ; then
    rm -rf ${wkdir}/${gene_name}_GReX/
    rm -f ${wkdir}/ABCA7_grex.geno
fi

exit

