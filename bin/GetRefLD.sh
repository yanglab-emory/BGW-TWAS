#!/usr/bin/bash

geno_dir=$1
pheno=$2
filehead=$3
LDdir=$4
Score_dir=$5
LDwindow=$6

module load tabix

bfGWAS_SS_dir="/home/jluningham/Projects/bfGWAS_SS"
# geno_dir=/home/jyang/Collaborations/IrwinSAGE/BU_GWASs_CDSymptomDimensions/lddetect_Dosage_Files

echo ${SGE_TASK_ID}
echo ${filehead}

line=$(head -n $SGE_TASK_ID $filehead | tail -n1)

echo ${line}

##### Set up temperary directory
TMP_NAME=`/usr/bin/mktemp -u XXXXXXXX`
# For array jobs, include the SGE task ID (in addition to the job ID):
TEMPDIR="/scratch/${JOB_ID}_${SGE_TASK_ID}_${TMP_NAME}"

if [ -e $TEMPDIR ]; then
  echo "Error. temp dir already exists on `hostname`"
  exit
else
   mkdir $TEMPDIR
fi

##### run BFGWAS
cd ${TEMPDIR}

echo Run BFGWAS with genotype: ${geno_dir}/${line}.vcf.gz
echo And phenotype: ${pheno} 

if [ -f ${LDdir}/${line}.LDcorr.txt.gz ] ; then
	echo ${LDdir}/${line}.LDcorr.txt.gz exists! 
	LDwindow=1
fi

echo LDwindow is $LDwindow

### With input genotype file in dosage format
${bfGWAS_SS_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${pheno} -maf 0 -o ${line} -LDwindow ${LDwindow} -GTfield DS -saveSS -zipSS

### With input genotype file in VCF format
# ${bfGWAS_SS_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${pheno} -maf 0 -o ${line} -LDwindow ${LDwindow} -saveSS -zipSS

echo Run BFGWAS Successfully to generate summary statistics!

##### Copy summary stat back
if [ ! -f ${LDdir}/${line}.LDcorr.txt.gz ] ; then
	echo Copy log.txt and LDcorr.txt.gz back to $LDdir
	cp -f  ${TEMPDIR}/output/${line}.LDcorr.txt.gz ${LDdir}/
	cp -f  ${TEMPDIR}/output/${line}.log.txt ${LDdir}/
else 
	echo Copy log.txt back to $Score_dir
	cp -f  ${TEMPDIR}/output/${line}.log.txt ${Score_dir}/
fi

echo Copy score.txt back $Score_dir
cp -f  ${TEMPDIR}/output/${line}.score.txt.gz ${Score_dir}/


rm -rf ${TEMPDIR}

exit



