# Bayesian Genome-wide (BGW) TWAS Software Usage

Tool **BGW-TWAS** is developed for leveraging both **cis-** and **trans-** eQTL based on a Bayesian variable selection model to predict genetically regulated gene expression (**GReX**). The product of estimated **cis-** and **trans-** eQTL effect sizes and their corresponding posterior causal probabilities (PCP) of being an eQTL will be taken as variant weights to conduct **TWAS** –– gene-based association test. This tool is implemented through several command-line scripts described in this manual. Example shell commands are provided in `./Run_BGW.sh`.

Please cite our BGW-TWAS paper if you use the tool:
>[*Bayesian Genome-wide TWAS Method to Leverage both cis- and trans-eQTL Information through Summary Statistics.* 2020 AJHG.](https://www.cell.com/ajhg/pdfExtended/S0002-9297(20)30291-3)

Please contact **Jingjing Yang (<jingjing.yang@emory.edu>)** if there is any issue.


---

- [Software Installation](#software-installation)
- [Input Files](#input-files)
	- [1. Training Gene Expression File](#1-training-gene-expression-file)
	- [2. Training VCF Genotype Files](#2-training-vcf-genotype-files)
	- [3. Test Individual-level GWAS Data](#3-test-individual-level-gwas-data)
- [Example Usage](#example-usage)
	- [Step 1. Obtain Summary Statistics](#step-1-obtain-summary-statistics)
	- [Step 2. Prune Genome Segments](#step-2-prune-genome-segments)
	- [Step 3. Training Gene Expression Prediction Model](#step-3-training-gene-expression-prediction-model)
	- [Step 4. Predict Bayesian GReX](#step-4-predict-bayesian-grex)
	- [Step 5. Association Studies](#step-5-association-studies)

---

## Software Installation

### 1. Compile **Estep_mcmc**
* Install required C++ libraries C++ libraries **zlib**, **gsl**, **eigen3**, **lapack**, **atlas**, **blas** are used to develop this tool. Please install these libraries to your system and include the library path `-I[path to libraries]` accordingly in the C++ compilation command line in the `Makefile`.

* Compile C++ library `./libStatGen/libStatGen.a` under your system by using the following commands:

```
cd BGW-TWAS/libStatGen/;
make clean;
make
```

* Compile C++ source code for the executible file *./bin/Estep_mcmc* that will be used to run the Estep MCMC algorithm to estimate eQTL effect sizes and the posterior causal probabilities (PCP) to be an eQTL, by using the following commands under `BGW-TWAS/` directory:

```
cd BGW-TWAS/;
make clean ;
make
```

* Even though a compiled executible file *./bin/Estep_mcmc* from our cluster is provided on GITHUB, please still compile one for your own system.

### 2. Additional Requirements
* Tool [**tabix**](https://www.htslib.org/doc/tabix.html)
* R library [**data.table**](https://github.com/Rdatatable/data.table/wiki/Installation)
* [GIT LFS](https://git-lfs.com/) is required to download VCF genotype files correctly for running the provided example commands. These should be around 20-50MB per file.

## Input Files

Input files of **BGW-TWAS** tool are all tab-seperated text files with lines seperated by `\n`:

* Reference transcriptomic data include gene expression file, genotype files of training samples, and a list of filenames for training genotype files, which are required for training Bayesian gene expression prediction models to estimate **cis-** and **trans-** eQTL effect sizes and their corresponding posterior causal probabilities of being an eQTL.
* Sample IDs in the gene expression file have to be the same as in VCF files.
* Note that genotype files can be of any names, and `cis-` or `trans-` blocks will be automatically determined based on the gene region as in the gene expression file. Each genotype file is suggested to include around `10,000` SNPs or `2.5MB` region, which can be determined based on the LD blocks of the same ethnicity.
* Individual-level GWAS data of test samples include genotype VCF files for test samples, a list of filenames for test genotype VCF files, and a phenotype file for test samples.
* Please also install [**Git LFS**](https://git-lfs.com/) for donwloading example large genotype VCF files correctly. Training VCF files are around 18MB and test VCF files are around 50MB.
* Summary-level GWAS data include a text file of Z-score statistics by single variant GWAS test.


### 1. Training Gene Expression File

* A file containing gene expression levels of training samples as in `./Example/ExampleData/Gene_Exp_example.txt`, with one gene per row, and one sample per column starting from the 6th column, columns seperated by tabs `
\t`. The first five columns are required to be gene annotation information including chromosome number, starting position, ending position, Gene ID, and Gene Name or second set of Gene ID. 
* Note that the GeneName in the fifth column will be used to extract gene expression data, which has to be provided. The fourth column could be the same as the fifth column.
* New line sperater should be of `\n` instead of `^M`. Otherwise, the file will not be read correctly.
* The gene expression levels in this file should be the residuals of a linear regression model that regresses out other confounding variables such as age, sex, top 5 genotype PCs from the `log2(normalized_read_counts)` of raw RNAseq data.
* Example row of the first 6 columns of gene `ABCA7` is as follows:

| CHROM |	GeneStart |	GeneEnd |	    TargetID     | GeneName |	 ROS20275399 |
| ----- | --------- | ------- | ---------------- | -------- | ------------ |
|  19	  |  1040101	| 1065571 |	 ENSG00000064687 |	 ABCA7  | 0.6707739044 |


### 2. Training VCF Genotype Files

* **[VCF Genotype files](http://samtools.github.io/hts-specs/VCFv4.1.pdf)** are required for training gene expression prediction models. The VCF genotype files should be one per genome block (variants of the same chromosome should be in the same block), sorted by position, and then zipped by `bgzip` (file names are of `[filehead].vcf.gz`), e.g., `./Example/ExampleData/genotype_data_files/Cis_geno_block.vcf.gz`. 
	- All VCF genotype files should be put under the same parent directory, `${geno_dir}`, which should also be provided. For example, see `./Example/ExampleData/genotype_data_files/*_geno_block.vcf.gz`.
	* **BGW-TWAS** tool will determine **cis-** or **trans-** SNPs based on the gene start/end information provided in the gene expression file.
	* Genotype files are supposed to be segmented based on LD. Genome blocks are expected to be approximately independent with ~5K-10K SNPs per block.
	* The same set of segmented genotype files work for all genes.
	* Example segmentation information derived by [LDetect](https://bitbucket.org/nygcresearch/ldetect/src/master/) are provided in `./Example/ExampleData/genotype_data_files/lddetect_1KG_*_hg19.bed ` for EUR, AFR, and ASN populations.

* **List of file heads** of VCF genotype files as in `./Example/ExampleData/geno_block_filehead.txt` is required. Each row of the list file is the file head of the VCF file of one genome block as in `[filehead].vcf.gz`. Note that the VCF file extension suffix `.vcf.gz` should not be included.


### 3. Individual-level GWAS Data for Test Samples

* **Test [VCF Genotype files](http://samtools.github.io/hts-specs/VCFv4.1.pdf)** are required for predicting GReX values of test samples that have individual-level GWAS data. Different from the training VCF genotype files, these test genotype files should be one per chromosome containing `_CHR[chr_num]_` in their file name, sorted by position, zipped by `bgzip`, and indexed by `tabix`. For example, see `./Example/ExampleData/genotype_data_files/Test_CHR*_geno.vcf.gz`.

* **Test phenotype file** contains two columns without header: the first column is a list of sample IDs that have genotype data in the test VCF files, and the second column is phenotype values. Second column will not be used for predicting GReX values, thus random values from `N(0, 1)` can also be used here without known phenotype values. For example, see `./Example/ExampleData/Test_pheno.txt`


## Example Usage

### Set up Tool Directories and Input Arguments:

```
## Variables for Step 1.

### Please update tool directory and use different working and LD file directories from below
BGW_dir=~/GIT/BGW-TWAS # tool directory
wkdir=${BGW_dir}/Example/ExampleWorkDir
LDdir=${BGW_dir}/Example/ExampleData/LDdir # shared by all genes
Score_dir=${BGW_dir}/Example/ExampleWorkDir/ABCA7_scores

### The following variables can stay the same for the testing purpose
GeneExpFile=${BGW_dir}/Example/ExampleData/Gene_Exp_example.txt
gene_name=ABCA7
geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
Genome_Seg_Filehead=${BGW_dir}/Example/ExampleData/geno_block_filehead.txt
GTfield=GT # specify genotype field "DS" or the corresponding field name in the VCF file for dosage genotype
num_cores=2 # number of cores to be used

## Varables for Step 2
p_thresh=0.001 # p-value threshold
max_blocks=50 # maximum blocks

## Variables for Step 3
N=499 # Training sample size
hfile=${BGW_dir}/Example/hypval.txt
PCP_thresh=0.0001

## Variables for Step 4
BGW_weight=${wkdir}/${gene_name}_BGW_eQTL_weights.txt
test_geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
test_geno_filehead=${BGW_dir}/Example/ExampleData/test_geno_filehead.txt
test_pheno=${BGW_dir}/Example/ExampleData/Test_pheno.txt
GTfield_test=GT # or DS if doseage data to be used
```

#### Remarks
* It is important to include the complete directory starting with `/home/` for all input files.
* Please set up working directory `${wkdir}` and LD file directory `${LDdir}` differently from the above examples, so that you may examine how your test outputs  be dimightfferent from the provided example outputs.
* `${LDdir}` stays the same for all genes.
* Please do not directly run the `./Run_BGW.sh` script. Instead, test the bash commands of the following steps one by one.

### Step 1. Obtain Summary Statistics
Shell script `Step1_get_sum_stat.sh` will generate single variant eQTL summary statistics (aka Score Statistics) in required formats.

#### Input arguments
- `--BGW_dir` : Specify the directory of BGW-TWAS tool
- `--wkdir` : Specify a working directory with writing access
- `-GeneExpFile` : Gene expression file directory
- `--gene_name` : Gene name that should be the same used in `GeneExpFile`
- `--geno_dir` : Directory of all genotype files
- `--LDdir` : Directory of all LD files with writing access
- `--Genome_Seg_Filehead` : Genome segmentation file
- `--GTfield` : Specify the genotype format in the training vcf file that should be used. Default `GT` for assayed genotype. Alternative value `DS` for the inputed dosage genotype field by [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!), or other genotype _FIELD_ name as used in the VCF file.
- `--num_cores` : Number of parallele sessions, default `1`.
- `--clean_output` : Whether to delete all intermediate outcomes, taking the input value of `1` for deleting or `0` keeping all intermediate files for testing purpose.

#### Example command:
```
# Change "--clean_output 0 " if you want to see intermediate files for debug purpose
${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} \
--Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores} --clean_output 1
```

* Rare variants with MAF <0.5% will be excluded from analysis by default.

#### Output files
* This shell script will create a directory called `${gene}_scores/` under the specified working directory `${wkdir}/`, where results for each genome segment will be stored.

* This script will also generate LD files under `${LDdir}/` for all genome blocks, if the LD files do not exist under `${LDdir}/`. Since LD files will be the same per genome block, these files will only be generated once and used for training prediction models for all gene expression traits.

* These summary statistics files and LD files will be used for implementing the *EM-MCMC* algorithm to fit the Bayesian model in Step 3.

* Gene expression phenotype file `${wkdir}/${gene_name}_exp_trait.txt` and variance file `${wkdir}/${gene_name}_exp_var.txt` will be used in Step 3.


### Step 2. Prune Genome Segments
Step 2 selects a subset of genome blocks (up to `${max_blocks}`) for joint model training by BGW-TWAS, where cis blocks are always selected. Trans blocks with minimum single-variant eQTL `p-value < ${p_thresh}` will first be ranked by the smallest p-value within block (from the smallest to the largest), where top ranked trans blocks will be selected up to `${max_blocks}`.


#### Input arguments
- `--wkdir` : Specify a working directory
- `--gene_name` : Gene name that should be the same used in `GeneExpFile`
- `-GeneExpFile` : Gene expression file directory
- `--Genome_Seg_Filehead` : Directory of the file containing a list of fileheads of segmented genotype files
- `--Score_dir` : Specify summary score statistic file directory, default `${wkdir}/${gene_name}_scores`
- `--p_thresh` : Specify p-value threshold for pruning, default `1e-5`
- `--max_blocks` : Specify the maximum number of genome blocks for jointly training gene expression prediction models, default `50`

#### Example command:
```
Score_dir=${wkdir}/${gene_name}_scores

# Change "--clean_output 0 " if you want to see intermediate files for debug purpose
${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--Score_dir ${Score_dir} --p_thresh ${p_thresh} --max_blocks ${max_blocks} \
--clean_output 1
```


#### Output files
* A list of filehead of selected genome blocks and the corresponding minimum within-block single eQTL analysis p-values are given (one genome block per row) in the file of `${wkdir}/${gene_name}_select_filehead.txt`, which will be used in step 3 (model training).

### Step 3. Training Gene Expression Prediction Model
Step 3 will use summary statistics generated from Step 1 for genome blocks pruned from Step 2 to carry out Bayesian variable selection regression with EM-MCMC algorithm, with the goal of fitting a gene expression prediction model. Bayesian posterior causal probability (PCP) to be an eQTL and effect sizes for SNPs with `PCP>${PCP_thresh}` will be saved for predicting GReX.

#### Input arguments
- `--BGW_dir` : Directory of BGW-TWAS tool
- `--wkdir` : Working directory with writting access
- `-GeneExpFile` : Gene expression file directory
- `--gene_name` : Study gene name that should be the same used in `GeneExpFile`
- `--LDdir` : Directory of all LD files with writting access
- `--Score_dir` : Specify summary score statistic file directory, default `${wkdir}/${gene_name}_scores`
- `--select_filehead` : File with selected genome block fileheads, default `${wkdir}/${gene_name}_select_filehead.txt`
- `--N` : Number of sample size used in Step 1 to generate eQTL summary statistics
- `--hfile` : Hyper parameter file as in `./Example/hypval.txt`, specifing the prior causal probability (_pi_) and effect size variance (_sigma2_) for cis (row 1) and trans (row 2) eQTL
- `--em` : Number of EM iterations. Default `3`.
- `--burnin` : Number of burnin MCMC iterations. Default `10000`.
- `--Nmcmc` : Number of MCMC iterations. Default `10000`.
- `--PCP_thresh` : PCP threshold for selecting eQTL with `PCP > ${PCP_thresh}` for predicting GReX and association study. Default `0.0001`.
- `--num_cores` : Specify the number of parallele sessions, default `1`.
- `--clean_output` : Whether to delete all intermediate outcomes, taking the input value of `1` for deleting or `0` keeping all intermediate files for testing purpose.

#### Example commands:
```
select_filehead=${wkdir}/${gene_name}_select_filehead.txt

# Change "--clean_output 0 " if you want to see intermediate files for debug purpose
${BGW_dir}/bin/Step3_EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --LDdir ${LDdir} \
--Score_dir ${Score_dir} --select_filehead ${select_filehead} \
--N ${N} --hfile ${hfile} \
--em 3 --burnin 10000 --Nmcmc 10000 \
--PCP_thresh ${PCP_thresh} --num_cores ${num_cores} \
--clean_output 1
```

* Intermediate output will be deleted unless with input argument `--clean_output 0`. Keeping intermediate outputs is recommended only for testing purpose.

#### Output files
* Bayesian estimates of eQTL PCP and effect sizes from the final EM-MCMC iteration will saved in the output `${gene_name}_BGW_eQTL_weights.txt` file under specified `${wkdir}` that lists all SNPs with `PCP>${PCP_thresh}`.

* The BGW weight file (`${gene_name}_BGW_eQTL_weights.txt`) will be used for predicting GReX values with individual-level GWAS data as in Step 4 or conducting gene-based association test with GWAS summary statistics. See instructions in Step 5 for TWAS procedure.

[//]: <> (Numerous arguments can be used to modify the EM-MCMC algorithm in Step 3, but should be done with caution. These arguments are detailed in Yang et al. 2017 https://github.com/yjingj/bfGWAS/blob/master/bfGWAS_Manual.pdf)


### Step 4. Predict Bayesian GReX
Step 4 will use the BGW weight file (`${gene_name}_BGW_eQTL_weights.txt`) generated from Step 3 and provided individual-level GWAS data to predict GReX values for test samples. Both test genotype VCF files (saved per chromosome, each file name containing the pattern of `_CHR[chr_num]_`) and test phenotype file should be provided.

The product of eQTL PCP and effect size will give an expected eQTL effect size that will be used as SNP weight for estimating the GReX values in Step 4.

#### Input arguments
- `--BGW_dir` : Directory of BGW-TWAS tool
- `--wkdir` : Working directory with writting access
- `--gene_name` : Study gene name that should be the same used in `GeneExpFile`
- `--BGW_weight` : Directory of the BGW eQTL weight file
- `--test_geno_dir` : Directory of all test genotype VCF files
- `--test_geno_filehead` : List of fileheads for test genotype VCF files
- `--test_pheno` : Directory of the test phentoype file
- `--GTfield` : Specify the genotype field in the VCF file to be used. Default `GT` for assayed genotype. Alternative value `DS` for the inputed dosage genotype field by [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!).
- `--num_cores` : Specify the number of parallele sessions, default `1`.
- `--clean_output` : Whether to delete all intermediate outcomes, taking the input value of `1` for deleting or `0` keeping all intermediate files for testing purpose.


#### Example commands:

```
${BGW_dir}/bin/Step4_get_test_grex.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--BGW_weight ${BGW_weight} --test_geno_dir ${test_geno_dir} \
--test_geno_filehead ${test_geno_filehead} \
--GTfield ${GTfield_test} --test_pheno ${test_pheno} \
--num_cores ${num_cores}
```

* Intermediate output will be deleted unless with input argument `--clean_output 0`. Keeping intermediate outputs is recommended only for testing purpose.


#### Output
* Predicted GReX values are saved in `${wkdir}/${gene_name}_pred_grex.txt`, which can be used to calculate prediction R2 and test the association between _GReX_ and _Phenotype of Interest_ (i.e., TWAS).

* File `${wkdir}/${gene_name}\_sumPCP` contains the sum of posterior causal probabilities (PCP) of all analyzed SNPs, the sum of cis-SNPs, and the sum of trans-SNPs, which are the expected number of total eQTL, cis-eQTL, and trans-eQTL for the target gene.

### Step 5. Association Studies

* If using summary-level GWAS data for TWAS, Step 4 will not need to be implemented. One can use the weight files `${wkdir}/${gene_names}_BGW_eQTL_weights.txt` to select the test SNPs presenting in this weight file. 

* Note that the eQTL weights estimated BGW are using only centered gene expressions and genotype data but not standardized to have standard deviation `1`. For a valid TWAS using GWAS summary statistics, the S-PrediXcan test statistic with genotype covariance matrix in the denominator must be used, as the FUSION test statistic with genotype correlation matrix in the denominator assumes eQTL weights are derived from standardized gene expressions and genotype data.

* With summary-level GWAS Z-score statistics <img src="https://render.githubusercontent.com/render/math?math=Z_i, i=1, \cdots, m"> for single SNP tests, eQTL effect sizes (i.e., weights `w`) estimated by BGW-TWAS method for `m` test SNPs
	*  "S-PrediXcan" TWAS Z score statistic is given by
<img src="https://render.githubusercontent.com/render/math?math=\frac{\sum_{i=1}^m \left(w_i \sigma^2_i Z_i \right) } {\sqrt{\mathbf{w'\Sigma w}}}">, with reference LD **Covariance matrix <img src="https://render.githubusercontent.com/render/math?math=\mathbf{\Sigma}">**, and reference genotype variance (diagonal values of the LD covariance matrix) of SNP `i` <img src="https://render.githubusercontent.com/render/math?math=\sigma^2_i">.
	* "FUSION" TWAS Z score statistic is given by <img src="https://render.githubusercontent.com/render/math?math=\frac{\sum_{i=1}^m \left(w_i Z_i\right) } {\sqrt{\mathbf{w'Vw}}}">, with reference LD **Correlation matrix `V`**.
	* One can use "TABIX" tool to extract reference genotype data of the test SNPs from external reference VCF files, for the purpose of calculating the reference covariance or correlation matrix more easily.

* Or one can use our other [TIGAR](https://github.com/yanglab-emory/TIGAR) tool to obtain TWAS results, where both "S-PrediXcan" and "FUSION" TWAS Z score tests were implemented.



* If using individual-level GWAS data for TWAS, Step 4 will need to be implemented to obtain predicted GReX values for all test samples. Then a simple single variant association test between _GReX_ and _Phenotype of Interest_ will result in the TWAS results.


## Options to save disk storage

* By default setting, temporary directories such as `${wkdir}/${gene_name}_scores/output/`, `${wkdir}/${gene_name}_EM_MCMC/`, `${wkdir}/${gene_name}_GReX/`, and intermediate files `${wkdir}/${gene_name}_exp_trait.txt`, `${wkdir}/${gene_name}_exp_var.txt`, `${wkdir}/ABCA7_grex.geno` would be deleted. If not, you can delete those.

* All score statistics files under `${wkdir}/${gene_name}_scores/` contain all single variant eQTL analysis test summary statistics, which is optional to be either saved for other usage or deleted.

* The weight files `${wkdir}/${gene_names}_BGW_eQTL_weights.txt` are the only file needed for TWAS, which can be gzipped to save storage.

