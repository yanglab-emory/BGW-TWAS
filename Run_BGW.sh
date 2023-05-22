
################################################################
################################################################
# Step 0: Set up directories of the tool, input files, and working directory for implementing BGW-TWAS.
################################################################
################################################################

### Varables for Step 1
# Tool directory
BGW_dir=~/GIT/BGW-TWAS

# Gene expression file has the following gene information in the first 5 columns:
# CHR, GeneStart, GeneEnd, TargetID/GeneID_1, GeneName/GeneID_2
# And gene expression levels from column 6 with samples in columns and genes in rows.
GeneExpFile=${BGW_dir}/Example/ExampleData/Gene_Exp_example.txt

# Specify gene name/identifier as in the 5th column of the gene expression file.
gene_name=ABCA7

# Parent directory of genotype files (VCF)
geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
# File with fileheads of all VCF files as in ${geno_dir}/[filehead].vcf.gz of all genome blocks
Genome_Seg_Filehead=${BGW_dir}/Example/ExampleData/geno_block_filehead.txt
# Specify the genotype field "GT" (called genotype) or "DS" (imputation dosage) to be used from the VCF files
GTfield=GT

# Working directory to write intermediate files and output files, unique per gene
wkdir=${BGW_dir}/Example/ExampleWorkDir
# wkdir=/home/jyang/GIT/BGW-TWAS/Test

# Parent directory of all LD files
LDdir=${BGW_dir}/Example/ExampleData/LDdir
Score_dir=/home/jyang/GIT/BGW-TWAS/Example/ExampleWorkDir/ABCA7_scores

# Number of cores/parallele_jobs to be used/implemented
num_cores=2

### Varables for Step 2
p_thresh=0.001 # p-value threshold
max_blocks=50 # maximum blocks

### Variables for Step 3
N=499 # sample size
hfile=${BGW_dir}/Example/hypval.txt
PCP_thresh=0.0001

### Variables for Step 4
BGW_weight=${wkdir}/${gene_name}_BGW_eQTL_weights.txt
test_geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
test_geno_filehead=${BGW_dir}/Example/ExampleData/test_geno_filehead.txt
test_pheno=${BGW_dir}/Example/ExampleData/Test_pheno.txt
GTfield_test=GT #or DS

################################################################
################################################################
# Step 1: obtain summary statistics (i.e., Score Statistics)
# Run single-variant eQTL analyses on the gene expression and genotypes data for the same training samples
################################################################
################################################################

${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores} --clean_output 1

################################################################
################################################################
# Step 2: Prune blocks
# Select a subset of genome blocks (up to ${max_blocks}) for joint model training by BGW-TWAS
# Cis blocks are always selected
# Trans blocks with minimum single-variant eQTL p-value < ${p_thresh} will be ranked by the smallest p-value within block (from the smallest to the largest).
# Top ranked trans blocks will be selected up to ${max_blocks}.
################################################################
################################################################
Score_dir=${wkdir}/${gene_name}_scores

${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--Score_dir ${Score_dir} \
--p_thresh ${p_thresh} --max_blocks ${max_blocks}

################################################################
################################################################
# Step 3: Training BGW-TWAS/BVSR gene expression prediction model by EM-MCMC algorithm
################################################################
################################################################
select_filehead=${wkdir}/${gene_name}_select_filehead.txt

${BGW_dir}/bin/Step3_EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --LDdir ${LDdir} \
--select_filehead ${select_filehead} \
--N ${N} --hfile ${hfile} \
--em 3 --burnin 10000 --Nmcmc 10000 \
--PCP_thresh ${PCP_thresh} --num_cores ${num_cores}

################################################################
################################################################
# Step 4: Predict GReX for test samples with individual-level GWAS data
################################################################
################################################################

${BGW_dir}/bin/Step4_get_test_grex.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--BGW_weight ${BGW_weight} --test_geno_dir ${test_geno_dir} \
--test_geno_filehead ${test_geno_filehead} \
--GTfield ${GTfield} --test_pheno ${test_pheno} \
--num_cores ${num_cores}

# optional script for removing files
