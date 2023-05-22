/*
	Bayesian Functional GWAS with Summary Statistics --- MCMC (bfGWAS_SS:MCMC)
    Copyright (C) 2017  Jingjing Yang

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


#ifndef __SS_H__                
#define __SS_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <bitset>
#include <cstring>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_integration.h"

#include "gzstream.h"   
#include "param.h"
//#include "io.h"
#include "ReadVCF.h"
#include "mathfunc.h"
#include "lapack.h"
#include "lm.h"


using namespace std;

typedef short int int16;
typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;


class CALCSS{
public: 
    bool zipSS;
    long int LDwindow;

    size_t UnCompBufferSize;
    vector <size_t> CompBuffSizeVec;
    bool Compress_Flag;
    vector<pair<int, double> > UcharTable;

    string file_out;
    size_t ni_total, ni_test;   //number of individuals
    size_t ns_total, ns_test;   //number of snps
    size_t n_type;
    vector<bool> indicator_idv; //indicator for individuals (phenotypes),1 analyzed
    vector<bool> indicator_snp; //indicator for SNPs: 1 analyzed

    // Phenotype information
    double pheno_mean, pheno_var, trace_G;

    vector<SNPINFO> snpInfo; // SNP information
    vector<SNPPOS> snp_pos;
    vector<double> SNPmean;
    

    //functions
    void CopyFromParam (PARAM &cPar);

    void GetSS(uchar **X, gsl_vector *y, vector< vector<double> > &LD, vector<double> &beta, vector<double> &beta_SE, vector<double> &U_STAT, vector<double> &SQRT_V_STAT, vector<double> &pval, vector<pair<size_t, double> > &pos_ChisqTest, vector<double> &xtx_vec, vector<double> &snp_var_vec, vector<double> &ni_effect_vec);

    void WriteSS(const vector< vector<double> > &LD, const vector<double> &beta, const vector<double> &beta_SE, const vector<double> &U_STAT, const vector<double> &SQRT_V_STAT, const vector<double> &pval);

};

// convert cov_ij to r2_ij
double Conv_xtx2_r2(const double &xtx_ij, const vector<double> &xtx_vec, const size_t &i, const size_t &j); 

void getXty(const vector<double> &beta, const vector<double> &xtx, vector <double> &Xty);

void getPval(const vector<double> &beta, const vector<double> &beta_sd, vector <double> &pval, vector<pair<size_t, double> > &pos_ChisqTest);

double getXtX(const vector< vector<double> > &LD, const size_t &pos_i, const size_t &pos_j, const vector<double> &xtx_vec);

double CalcResVar(const gsl_vector * Xty_cond, const gsl_vector * beta_cond, const double &yty); 

void CalcBeta(const gsl_matrix *XtX_cond, const gsl_vector * Xty_cond, gsl_vector * beta_cond);


#endif


