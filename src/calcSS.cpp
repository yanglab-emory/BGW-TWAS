/*
	Bayesian Functional GWAS with Summary Statistics --- MCMC (BFGWAS_SS:MCMC)
    Copyright (C) 2018  Jingjing Yang

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


#include "calcSS.h"


void CALCSS::CopyFromParam (PARAM &cPar) 
{

    zipSS=cPar.zipSS;

    UnCompBufferSize = cPar.UnCompBufferSize;
    CompBuffSizeVec = cPar.CompBuffSizeVec;
    Compress_Flag = cPar.Compress_Flag;

    file_out=cPar.file_out;
        
    ni_total=cPar.ni_total;
    ns_total=cPar.ns_total;
    ni_test=cPar.ni_test;
    ns_test=cPar.ns_test;
    n_type = cPar.n_type;

    LDwindow=cPar.LDwindow;
    
    indicator_idv=cPar.indicator_idv;   
    indicator_snp=cPar.indicator_snp;

    SNPmean = cPar.SNPmean; 
    pheno_mean = cPar.pheno_mean;
    pheno_var = cPar.pheno_var;
    snp_pos = cPar.snp_pos;
    
    return;
}


//calculat summary statistics of score statistics and LD matrix
void CALCSS::GetSS(uchar **X, gsl_vector *y, vector< vector<double> > &LD, vector<double> &beta, vector<double> &beta_SE, vector<double> &U_STAT, vector<double> &SQRT_V_STAT, vector<double> &pval, vector<pair<size_t, double> > &pos_ChisqTest, vector<double> &xtx_vec, vector<double> &snp_var_vec, vector<double> &ni_effect_vec){

    cout << "\nStart calculating summary statistics ... \n";

    double yty;
    // Center y is centered by cPar.CopyPheno()
    gsl_blas_ddot(y, y, &yty); 
    pheno_var = yty / ((double)(ni_test-1)) ;
    cout << "ni_test = " << ni_test << endl;
    cout << "ns_test = " << ns_test << endl;
    cout << "pheno_var = " << pheno_var << "\n";

    //cout << "create UcharTable ...\n";
    CreateUcharTable(UcharTable);

    // define used variables 
    gsl_vector *xvec_i = gsl_vector_alloc(ni_test);
    gsl_vector *xvec_j = gsl_vector_alloc(ni_test);
    gsl_vector *xbeta_i = gsl_vector_alloc(ni_test);

    double xtx_ij, xty, xtx_i, beta_i, v2, chisq_i, beta_SE_i, r2;

    // cout << "calculate xtx by the order of chr/bp ... \n";
    xtx_vec.clear();
    snp_var_vec.clear();
    ni_effect_vec.clear();
    trace_G = 0;
    for (size_t i=0; i<ns_test; ++i) {
        //calculate xtx_i
        getGTgslVec(X, xvec_i, snp_pos[i].pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
        gsl_blas_ddot(xvec_i, xvec_i, &xtx_i);
        trace_G += xtx_i;
        xtx_vec.push_back( xtx_i );
        snp_var_vec.push_back (xtx_i / double(ni_test) );
        ni_effect_vec.push_back(ni_test);
    }

    // cout << "calculate beta, score statistics by the order of chr/bp ... \n";
    beta.clear();
    LD.clear(); 
    pval.clear();
    pos_ChisqTest.clear();
    beta_SE.clear();
    for (size_t i=0; i<ns_test; ++i) {
        //calculate xtx_i
        getGTgslVec(X, xvec_i, snp_pos[i].pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
        xtx_i = xtx_vec[i];

        //calculate effect-size
        if(xvec_i->size != y->size){cerr << "Genotype length dose not equal to phenotype length!\n Some samples in the genotype file may not have genotype data!\n Please check your phenotype and genotype input files!\n"; exit(-1);}
        gsl_blas_ddot(xvec_i, y, &xty);
        if(xtx_i > 0) beta_i = xty / xtx_i;
        else beta_i = 0.0;
        beta.push_back(beta_i); // effect size
        U_STAT.push_back(xty); // score statistic 
        v2 = pheno_var * xtx_i ;
        SQRT_V_STAT.push_back( sqrt(v2) ); // score statistic standard deviation
        chisq_i = ((double)ni_test)*(log(yty)-log(yty-xty*xty/xtx_i)); // LRT statistic
        // chisq_i = xty * xty / v2; // Score test statistic
        pval.push_back( gsl_cdf_chisq_Q (chisq_i, 1.0) ); // pvalue needed for BVSRM
        pos_ChisqTest.push_back( make_pair(i, chisq_i) ) ; // pos_ChisqTest needed for BVSRM

        gsl_vector_memcpy(xbeta_i, xvec_i);
        gsl_vector_scale(xbeta_i, -beta_i);
        gsl_vector_add(xbeta_i, y);
        gsl_blas_ddot(xbeta_i, xbeta_i, &beta_SE_i); // effect-size deviation
        if(xtx_i > 0) beta_SE_i = sqrt( beta_SE_i / ((double)ni_test * xtx_i) );
        else beta_SE_i = 0.0; 
        beta_SE.push_back(beta_SE_i);
        
        // saving X'X to LD
        LD.push_back(vector<double>()); // save correlation
        LD[i].push_back(snp_var_vec[i]); // save snp genotype variance

        if(i < (ns_test-1) ){
            //calculate xtx_ij 
            for(size_t j=(i+1); j < ns_test; ++j){
                if( (snp_pos[j].chr == snp_pos[i].chr) && (snp_pos[j].bp <= snp_pos[i].bp + LDwindow) )
                {
                    getGTgslVec(X, xvec_j, snp_pos[j].pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
                    gsl_blas_ddot(xvec_i, xvec_j, &xtx_ij);
                    r2 = Conv_xtx2_r2(xtx_ij, xtx_vec, i, j); 
                    LD[i].push_back( r2 ); // Correlation between x_i and x_j
                }
                else{break;}
            }
        }
    }

    gsl_vector_free(xvec_i);
    gsl_vector_free(xvec_j);
    gsl_vector_free(xbeta_i);

    return;
}

double Conv_xtx2_r2(const double &xtx_ij, const vector<double> &xtx_vec, const size_t &i, const size_t &j){
    // get correlation as r2
    double r2 = 0.0;
    if(xtx_ij != 0.0){
        if( (xtx_vec[i] > 0.0) && (xtx_vec[j] > 0.0) ){
            r2 = xtx_ij / sqrt( xtx_vec[i] * xtx_vec[j] );
        }
    }
    return r2;
}


void CALCSS::WriteSS(const vector< vector<double> > &LD, const vector<double> &beta, const vector<double> &beta_SE, const vector<double> &U_STAT, const vector<double> &SQRT_V_STAT, const vector<double> &pval)
{
    cout << "\nStart writing summary statistics ... \n";
    String fout = file_out.c_str();

    // output files matches RareMetalWorker outputs
    String cov_file_str = "./output/" + fout;
    String score_file_str = "./output/" + fout;

    IFILE cov_out=NULL;
    IFILE score_out=NULL;

    if(zipSS){
        cov_file_str +=".LDcorr.txt.gz";
        cov_out = ifopen(cov_file_str, "w", InputFile::BGZF);

        score_file_str += ".score.txt.gz";
        score_out = ifopen(score_file_str, "w", InputFile::BGZF);

        if(cov_out == NULL || score_out == NULL){
            perror("Fail to open LD or beta file!!! \n");
        }
    }else{
        cov_file_str +=".LDcorr.txt";
        cov_out = ifopen(cov_file_str, "w", InputFile::UNCOMPRESSED);

        score_file_str += ".score.txt";
        score_out = ifopen(score_file_str, "w", InputFile::UNCOMPRESSED);

        if(cov_out == NULL || score_out == NULL){
            perror("Fail to open LD or beta file!!! \n");
        }
    }

    // write an extra column saving xtx with centered genotypes
    ifprintf(score_out, "#CHROM\tPOS\tID\tREF\tALT\tN\tMAF\tHWE_PVALUE\tU_STAT\tSQRT_V_STAT\tEFFSIZE_BETA\tBETA_SE\tPVALUE\n");
    // assuming variants have unique CHR:POS 
    ifprintf(cov_out, "#CHROM\tPOS\tID\tREF\tALT\tN\tMAF\tCORR\n");
    
    //Write files by the order of chr/bp
    for(size_t i=0; i<ns_test; i++){

        // write score statistics
        ifprintf(score_out, "%s\t%ld\t%s\t%s\t%s\t%u\t%g\t%s\t%g\t%g\t%g\t%g\t%g\n", snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].rs.c_str(), snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str(), ni_test, snp_pos[i].maf, "NA", U_STAT[i], SQRT_V_STAT[i], beta[i], beta_SE[i], pval[i]);

        // write banded covariance matrix: chr pos ref alt
        ifprintf(cov_out, "%s\t%ld\t%s\t%s\t%s\t%u\t%g\t", snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].rs.c_str(), snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str(), ni_test, snp_pos[i].maf);

      //  for(size_t j=0; j<LD[i].size(); j++){
      //      ifprintf(cov_out, "%ld,", snp_pos[i+j].bp);
      //  }
      //  ifprintf(cov_out, "\t");

        for(size_t j=0; j<LD[i].size(); j++){
            ifprintf(cov_out, "%g,", LD[i][j]);
        }
        ifprintf(cov_out, "\n");

    }

    ifclose(cov_out);
    ifclose(score_out);

    // tabix zipped files
    String cmd;
    int sys_status=1;

    if(zipSS){
        printf("Tabixing .LDcorr.txt.gz files ... \n");
        cmd = String("tabix -c \"#\" -s 1 -b 2 -e 2 -f ") + cov_file_str;
        sys_status = system(cmd.c_str());
        if ( sys_status == 0 ) {
            printf( "LD correlation output %s has been tabixed\n", cov_file_str.c_str() );
        }
        else {
            printf("Unable to tabix %s\n", cov_file_str.c_str());
        }

        printf("Tabixing .score.txt.gz files ... \n");
        cmd = String("tabix -c \"#\" -s 1 -b 2 -e 2 -f ") + score_file_str;
        sys_status = system(cmd.c_str());
        if ( sys_status == 0 ) {
            printf( "Score statistics output %s has been tabixed\n", score_file_str.c_str() );
        }
        else {
            printf("Unable to tabix %s\n", score_file_str.c_str());
        }
    }

    return;
}

void getXty(const vector<double> &beta, const vector<double> &xtx, vector <double> &Xty)
{
	//n is the sample size
    cout << "Calculate Xty ... \n";
    Xty.clear();
    for(size_t i=0; i<beta.size(); i++){
        Xty.push_back( beta[i] * xtx[i] );
    }
    return;
}

void getPval(const vector<double> &beta, const vector<double> &beta_sd, vector <double> &pval, vector<pair<size_t, double> > &pos_ChisqTest)
{
    cout << "Calculate pval ... \n";
    pval.clear();
    pos_ChisqTest.clear();
    double pval_i, chisq_i;

    for(size_t i=0; i<beta.size(); i++){
        chisq_i = pow(beta[i] / beta_sd[i], 2);
        pos_ChisqTest.push_back( make_pair(i, chisq_i) );

        pval_i = gsl_cdf_chisq_Q (chisq_i, 1);
        pval.push_back(pval_i);
    }
    return;
}

// LD has been set up for analyzed variants
// xtx_vec has been setup based on snp variance and effective sample size
double getXtX(const vector< vector<double> > &LD, const size_t &pos_i, const size_t &pos_j, const vector<double> &xtx_vec)
{
    double xtx_ij = 0.0;

    if(pos_i == pos_j){
        xtx_ij = xtx_vec[pos_i];
    }
    else 
    {
        if( (pos_j - pos_i) > 0 && (pos_j - pos_i) < LD[pos_i].size()  ) 
            {
                xtx_ij = LD[pos_i][pos_j - pos_i] * sqrt(xtx_vec[pos_i] * xtx_vec[pos_j]);   
            }     
        else if( (pos_i - pos_j) > 0 && (pos_i - pos_j) < LD[pos_j].size() ) 
            {
                xtx_ij = LD[pos_j][pos_i - pos_j] * sqrt(xtx_vec[pos_i] * xtx_vec[pos_j]);
            }
    }

    return xtx_ij;
}



double CalcResVar(const gsl_vector * Xty_cond, const gsl_vector * beta_cond, const double &yty)
{
    double rtr, xtyb;

    gsl_blas_ddot(Xty_cond, beta_cond, &xtyb);

    //cout << "Regression R2 in calcResVar = " << xtyb / yty << endl;

    rtr = yty - xtyb ;

    if(rtr <= 0){
        cout << "Regression R2 in calcResVar = " << xtyb / yty << endl;
        perror("Nonpositive residual variance!\n");  
    }

    return rtr;
}


void CalcBeta(const gsl_matrix *XtX_cond, const gsl_vector * Xty_cond, gsl_vector * beta_cond)
{
    size_t s_size = Xty_cond->size;

    gsl_matrix *XtXinv = gsl_matrix_alloc(s_size, s_size);

    gsl_matrix_memcpy(XtXinv, XtX_cond);

    LapackSolve(XtXinv, Xty_cond, beta_cond);

    gsl_matrix_free(XtXinv);

    return ;
}












