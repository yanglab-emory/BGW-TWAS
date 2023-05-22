/*
 Bayesian Functional GWAS --- MCMC (bfGWAS:MCMC)
Copyright (C) 2016  Jingjing Yang

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __PARAM_H__                
#define __PARAM_H__

#include <vector>
#include <map>
#include <set>
#include <iostream>     // std::cout, std::endl
#include <iomanip>
#include <limits>


#include "VcfFileReader.h"
#include "MemoryAllocators.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cstring>
#include <sys/stat.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"

#include "lapack.h"
#include "gzstream.h"
#include "mathfunc.h"
#include "ReadVCF.h"


using namespace std;

typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;


class SNPINFO {
public:
  string chr;
  string rs_number;
  double cM;
  long int base_position;
  string a_minor;
  string a_major;
  int n_miss;
  double missingness;
  double maf;
  
  vector<bool> indicator_func;
  vector<double> weight;
  double weight_i;
  
  string key;
  
  void printMarker();
};

struct genMarker
{
  string rs;
  string chr;
  long int bp;
  string Ref;
  string Alt;
  
  void iniRecord(VcfRecord& record);
  void printMarker();
  
};

//JY
struct SNPPOS{
  size_t pos;
  string rs;
  string chr;
  long int bp;
  string a_minor;
  string a_major;
  double maf;
  vector<bool> indicator_func;
  string key;
  //vector<double> weight;
  //double weight_i;
  
  void printMarker();
};

struct AnnoINFO {
public:
  string chr;
  string rs_number;
  long int base_position;
  vector<bool> indicator_func;
  vector<double> weight;
  double weight_i;
};


//JY
void printSNPInfo(vector<SNPPOS> &snp_pos, int numMarker);
void CalcWeight(const vector<bool> &indicator_func, vector<double> &weight, const double weight_i);

//results for lmm
class SUMSTAT {
public:
  double beta;				//REML estimator for beta
  double se;				//SE for beta  
  double lambda_remle;		//REML estimator for lambda
  double lambda_mle;		//MLE estimator for lambda
  double p_wald;			//p value from a Wald test
  double p_lrt;				//p value from a likelihood ratio test
  double p_score;			//p value from a score test
  
  double score_u; // for meta-analysis
  double score_v;
};

//results for mvlmm
class MPHSUMSTAT {
public:
  vector<double> v_beta;	//REML estimator for beta
  double p_wald;			//p value from a Wald test
  double p_lrt;				//p value from a likelihood ratio test
  double p_score;			//p value from a score test
  vector<double> v_Vg;	//estimator for Vg, right half
  vector<double> v_Ve;	//estimator for Ve, right half
  vector<double> v_Vbeta;	//estimator for Vbeta, right half
};


//hyper-parameters for bslmm
class HYPBSLMM {
public:
  double h;
  double rho;
  double logp;
  
  vector<double> theta;
  vector<double> log_theta; // log(theta) for each function type
  vector<double> subvar; // variance for each function type
  vector<double> rho_vec;
  vector<size_t> m_gamma; // # of selected SNP of each function type
  
  double rv;
  double sigma_b2; //
  double pve;
  double pge;
  size_t n_gamma;
  
};


class PARAM {
public:
  //multiple functionrelated parameters
  size_t n_type;
  vector<size_t> mFunc; // # of variants of each variant type
  double e; //hyper parameter in the prior gamma distribution
  double vscale;
  map<string, int> mapFunc2Code;
  int iniType;
  bool FIXHYP;
  
  bool calc_K, saveGeno, saveSS, zipSS, inputSS, refLD, scaleN, printLD, use_xtx_LD, final_EM;
  long int LDwindow;
  
  string iniSNPfile;
  string hypfile;
  double rv, pheno_mean, pheno_var;
  
  // IO related parameters
  size_t UnCompBufferSize;
  vector <size_t> CompBuffSizeVec;
  bool Compress_Flag;
  map<string, size_t>  PhenoID2Ind; // map sampleID to index in y (for reorder pheno with respect to geno)
  map<string, size_t>  PhenoID2Pos; // map sampleID to index in raw pheno (for obtaining pheno_index to be used with indicator_idv)
  map<string, size_t> GenoSampleID2Ind; // index sampleIDs with both pheno and geno data
  vector<size_t> SampleVcfPos;
  vector<string> VcfSampleID; // size=total sample #
  vector<string> VcfSampleID_test; // sample id for ni_test in order of X[i][j]
  vector<string> InputSampleID; //size = ni_total
  vector<pair<int, double> > UcharTable;
  vector<double> SNPmean;
  
  // JML
  // Annotation-related parameters
  string target_chr;
  long int start_pos; 
  long int end_pos;
  long int window_size;
  
  
  // SS related parameters
  string file_score;
  string file_cov;
  vector< vector<double> > LD_ref, LD;
  vector<double> U_STAT, SQRT_V_STAT, pval_vec, xtx_vec, snp_var_vec;
  vector<double>  mbeta, mbeta_SE, ni_effect_vec;
  vector<pair<size_t, double> >  pos_ChisqTest;
  
  bool mode_silence;
  int a_mode;				//analysis mode, 1/2/3/4 for Frequentist tests
  int k_mode;				//kinship read mode: 1: n by n matrix, 2: id/id/k_value; 
  size_t d_pace;		
  
  string file_func_code; //coded all unique variant function types
  string file_sample; // file containing analyzed sample IDs
  string file_vcf;
  string GTfield;
  
  string file_bfile;
  string file_geno;
  string file_pheno;
  string file_anno;		
  string file_kin;
  string file_out;
  string file_snps;		//file containing analyzed snps or genes
  string file_grex; // file of genotypes for predicting GREX in final repetition. 
  
  // run prediction in second command
  // second class for reading in param estimates from training 
  // so it can match up training scores and prediction data
  // 10/7/19
  
  // QC related parameters	
  double miss_level;
  double maf_level;	
  double hwe_level;
  double r2_level;
  
  // BVSRM MCMC related parameters
  size_t win;
  size_t nadd_accept, ndel_accept, nswitch_accept, nother_accept;
  size_t nadd, ndel, nswitch, nother;
  
  double h_min, h_max, h_scale;			//priors for h
  double rho_min, rho_max, rho_scale;		//priors for rho
  double logp_min, logp_max, logp_scale;		//priors for log(pi)
  size_t s_min, s_max;			//minimum and maximum number of gammas
  size_t w_step;					//number of warm up/burn in iterations
  size_t s_step;					//number of sampling iterations
  size_t n_accept;				//number of acceptance
  size_t n_mh;					//number of MH steps within each iteration
  size_t region_pip;              //number of MCMC with SNPs>0
  long int randseed;
  double trace_G;
  
  HYPBSLMM cHyp_initial;
  
  // JML: VB-related parameters
  
  size_t max_iter; //max VB iterations for convergence
  double convergence;
  // delta (convergence) level? 
  
  // is this the only new parameter needed? 
  
  
  // Summary statistics
  bool error;
  size_t ni_total, ni_test, ni_cvt;	//number of individuals
  size_t np_obs, np_miss;		//number of observed and missing phenotypes
  
  size_t ns_total, ns_test;	//number of snps
  size_t n_cvt;
  
  size_t ng_total, ng_test;	//number of genes
  size_t ni_control, ni_case;	//number of controls and number of cases
  double time_total;		//record total time
  double time_G;			//time spent on reading files the second time and calculate K
  double time_opt;		//time spent on optimization iterations/or mcmc
  double time_Omega;		//time spent on calculating Omega
  double time_hyp;		//time spent on sampling hyper-parameters, in PMM
  double time_Proposal;  //time spend on constructing the proposal distribution (i.e. the initial lmm or lm analysis)
  
  // Data
  
  vector<double>  pheno;	//a vector record all phenotypes, NA replaced with -9
  
  vector<bool>  indicator_pheno;			//a matrix record when a phenotype is missing for an individual; 0 missing, 1 available
  
  vector<bool> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis
  
  vector<bool> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
  
  vector<bool> indicator_idv_grex;       //indicator for individuals in GREX prediction
  
  vector<bool> indicator_snp_grex;       //sequence indicator for SNPs: 0 not in GREX set; 1 available for prediction

  
  map<string, size_t> mapScoreKey2Pos; //map snpkey to the position in score.txt
  map<string, size_t> mapLDKey2Pos; // map snpkey to the position in LDR2.txt
  map<int, string> mapCode2Func; // map unique code to a unique function type
  
  map<string, int> mapID2num;		//map chr:pos:ref:alt to position, from 0 to n-1
  map<string, string> mapRS2chr;		//map rs# to chromosome location
  map<string, long int> mapRS2bp;		//map rs# to base position
  map<string, double> mapRS2cM;		//map rs# to cM
  map<string, double> mapRS2est;			//map rs# to parameters

  // map model-building SNP to GREx SNP? 
  // basically, what happens if we used a SNP in the bvsr, it has > 0 pi, but it is not in the GREx sample? 

  
  vector<SNPINFO> snpInfo;		//record SNP information
  vector<SNPPOS> snp_pos;
  set<string> setSnps;			//a set of snps for analysis

  set<string> setPredSnps; //SNPs in the GREX prediction sample 
  
  //constructor
  PARAM();
  
  //functions
  void ReadFiles ();		
  void CheckParam (); 
  void CheckData ();	
  void PrintSummary ();
  void ReadGenotypes (uchar **X, gsl_matrix *K);
  void WriteGenotypes(uchar **X);
  
  // summary statistics;
  void ReadSS ();
  void UpdateScore();
  
  //void CheckCvt ();
  void ProcessPheno();
  void UpdatePheno();
  void ReorderPheno(gsl_vector *y);
  
  void CopyPheno (gsl_vector *y);
  void CalcKin (gsl_matrix *matrix_kin);
  
  // get LD from LDref
  void Convert_LD();
  
};


void CreateSnpPosVec(vector<SNPPOS> &snp_pos, vector<SNPINFO> &snpInfo, const vector<bool> &indicator_snp);

bool comp_snp(const SNPPOS& lhs, const SNPPOS& rhs);

// summary stat
vector<string> split(const string& str, const string& delim);

void SwapKey(string &key);


// get r2 from LDref
double getR2_ij(const vector< vector<double> > &LD, size_t idx_i, size_t idx_j, const bool &swap_i, const bool &swap_j);


#endif

