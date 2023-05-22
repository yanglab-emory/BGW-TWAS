/*
	Bayesian Functional GWAS --- MCMC (BFGWAS:MCMC)
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


#include "param.h"
#include "bvsrm.h"
#include "io.h"

// define to_string function : convert to string
template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}


void genMarker::iniRecord(VcfRecord& record){

    rs = record.getIDStr();
    chr = record.getChromStr();
    bp = record.get1BasedPosition();
    Ref = record.getRefStr();
    Alt = record.getAltStr();

    return;
}

void genMarker::printMarker(){

    std::cout << "ID : " << rs <<"; ";
    std::cout << "chr : " << chr <<"; ";
    std::cout << "bp : " << bp <<"; ";
    std::cout << "Ref : " << Ref <<"; ";
    std::cout << "Alt : " << Alt <<"\n";

    return;
}

void SNPPOS::printMarker(){
    std::cout << "position : " << pos << "; ";
    std::cout << "ID : " << rs <<"; ";
    std::cout << "chr : " << chr <<"; ";
    std::cout << "bp : " << bp <<"; ";
    std::cout << "minor allel : " << a_minor <<"; ";
    std::cout << "major allel : " << a_major <<"; \n";
}

void SNPINFO::printMarker(){

    std::cout << "ID : " << rs_number <<"; ";
    std::cout << "chr : " << chr <<"; ";
    std::cout << "bp : " << base_position <<"; ";
    std::cout << "Ref : " << a_major <<"; ";
    std::cout << "Alt : " << a_minor <<"\n";
    std::cout << "missingness = " << missingness<< "; maf = " << maf << "\n";
}

void printSNPInfo(vector<SNPPOS> &snp_pos, int numMarker)
{
    for (int i=0; i<numMarker; i++) {
        snp_pos[i].printMarker();
    }
}

void CalcWeight(const vector<bool> &indicator_func, vector<double> &weight, const double weight_i)
{
    weight.clear();
    for (size_t i=0; i < indicator_func.size(); i++) {
            if (indicator_func[i]) weight.push_back(weight_i);
            else weight.push_back(0.0);
        }
    if (weight.size() != indicator_func.size()) {
        cout << "Error weight size.\n";
    }
}


PARAM::PARAM(void):
vscale(0.0), iniType(3), calc_K(0), saveGeno(0), saveSS(0), zipSS(0),
inputSS(0), refLD(0), scaleN(0), printLD(0), use_xtx_LD(0), LDwindow(1000000), rv(0.0), Compress_Flag(0),
mode_silence (false), a_mode (0), k_mode(1), d_pace (100000),
GTfield("GT"), file_out("result"), final_EM(0),
miss_level(0.05), maf_level(0.005), hwe_level(0), r2_level(0.001),
win(100),nadd_accept(0), ndel_accept(0), nswitch_accept(0),
nother_accept(0), nadd(0), ndel(0),
nswitch(0), nother(0),
h_min(-1), h_max(1.0), h_scale(-1),
rho_min(1.0), rho_max(1.0),	rho_scale(-1),
logp_min(0.0), logp_max(0.0), logp_scale(-1),
s_min(0), s_max(10), max_iter(1000), convergence(1e-7),
w_step(50000),	s_step(500000), n_accept(0),
n_mh(10), randseed(2016), error(false), ni_test(0),
time_total(0.0), time_G(0.0), time_Omega(0.0),
target_chr("12"), start_pos(241229212), end_pos(241315842), window_size(1000000)
{}


bool comp_snp(const SNPPOS& lhs, const SNPPOS& rhs){
    return (lhs.chr.compare(rhs.chr) < 0) || ((lhs.chr.compare(rhs.chr) == 0) && (lhs.bp < rhs.bp));
}


//read individual level data files
//obtain ns_total, ng_total, ns_test, ni_test, n_type
void PARAM::ReadFiles (void)
{
	cout<<"\nStart reading individual-level data ...\n";
	string file_str;

	// read the set of to be analyzed SNP ids if given file_snps
	if (!file_snps.empty()) {
		if (ReadFile_snps (file_snps, setSnps)==false) {error=true;}
	} else {
		setSnps.clear();
	}

		//read bed/bim/fam files of plink format
		if (!file_bfile.empty()) {
			cout << "\nStart reading plink bim/fam files ...\n";
			file_str=file_bfile+".bim";
			if (ReadFile_bim (file_str, snpInfo, mapID2num)==false) {error=true;}

			file_str=file_bfile+".fam";
			if (ReadFile_fam (file_str, indicator_idv, pheno, InputSampleID, ni_total)==false) {error=true;}
			if(!file_sample.empty()){
				// revise indicator_idv based on analyzed samples
				readFile_sample (file_sample, InputSampleID, indicator_idv);
			}
			// set up VcfSampleID for fam file
			cout << "Set up VcfSampleID as analyzed sample IDs\n";
			VcfSampleID.clear();
			for(size_t i=0; i < indicator_idv.size(); i++){
				if(indicator_idv[i]){
					VcfSampleID.push_back(InputSampleID[i]);
				}
			}

			// obtain ni_test, ni_total, PhenoID2Ind, PhenoID2Pos before reading genotypes
			ProcessPheno();

			file_str=file_bfile+".bed";
			cout << "First time reading Plink bed file: \n";
			if (ReadFile_bed (file_str, setSnps, indicator_idv, indicator_snp, snpInfo, PhenoID2Pos, ni_test, ni_total, maf_level, miss_level, hwe_level, ns_test, ns_total)==false) {error=true;}
			//cout << "First time reading Plink bed file success! \n";
	    }else{
	    	if (!file_pheno.empty()){
	    		cout << "\nStart reading pheno file ...\n";
	        	if (ReadFile_pheno (file_pheno, indicator_idv, pheno, InputSampleID, ni_total)==false)
	            	{error=true;}
	            if(!file_sample.empty()){
					// revise indicator_idv based on analyzed samples
					readFile_sample (file_sample, InputSampleID, indicator_idv);
				}
	        	ProcessPheno();
	        	// obtain ni_test, ni_total, PhenoID2Ind, PhenoID2Pos before reading genotypes
	    	}else{
	    		cout << "No phenotype input file, extracting sample information from the vcf/genotype files.\n";
	    		if (!file_vcf.empty()) {
	        		getIDVvcf(file_vcf, indicator_idv,  ni_total, GTfield);
	      		}else if (!file_geno.empty()) {
					getIDVgeno(file_geno, indicator_idv, ni_total) ;
		  		}else{
		  			cerr << "Unable to get sample information!" << endl;
		  			exit(-1);
		  		}
	    	}

	      	//read vcf file for genotypes
	      	if (!file_vcf.empty()) {
	        	cout << "\nStart reading vcf file first time ...\n";
	        	indicator_snp.clear();
	        	snpInfo.clear();
	        	if (ReadFile_vcf(file_vcf, setSnps, indicator_idv, indicator_snp, maf_level, miss_level, hwe_level, snpInfo, ns_test, ns_total, ni_test, GTfield, PhenoID2Pos, VcfSampleID, SampleVcfPos, mapID2num, GenoSampleID2Ind) == false )
	            	{error=true;}
	      	}else if (!file_geno.empty()) {
		  		//read genotype file
		  		cout << "\nStart reading dosage file first time ...\n";
				if (ReadFile_geno (file_geno, setSnps, indicator_idv, indicator_snp, PhenoID2Pos, snpInfo, VcfSampleID, SampleVcfPos, maf_level, miss_level, hwe_level, ns_test, ns_total, ni_test, ni_total, mapID2num, GenoSampleID2Ind)==false) {error=true;}
		  	}

		  	// update ni_test, ni_total, PhenoID2Ind based on samples that also have geno data
			UpdatePheno();
		}

	    if ( (!file_anno.empty()) && (!file_func_code.empty()) ) {
	    	cout << "\nStart reading annotation files ...\n";
	    	//cout << file_anno << " \nwith code file " << file_func_code << "\n";
	        if (ReadFile_anno (file_anno, file_func_code, mapFunc2Code, indicator_snp, snpInfo, n_type, mFunc, mapID2num)==false) {error=true;}
	    }
	    else{
            cout << "\nAnnotation assignment by position!\n";
            if (setANNOCode (start_pos, end_pos, target_chr, window_size, indicator_snp, snpInfo, n_type, mFunc)==false) {error=true;}
	    }
    // For the future: clean this up by allowing empty annotation if target chr, position, and window_size are not provided!
   // else{
    //    cout << "\nEmpty annotation file, no target position, all variants are treated as of one category!\n";
     //   if (Empty_anno (indicator_snp, snpInfo, n_type, mFunc)==false) {error=true;}

   // }

	    // need Snp info

// jy edit on 11/2020
	   // if (EM_step=final_EM)
	   // if (!file_grex.empty()) //

	return;
}


void PARAM::CheckParam (void)
{
	struct stat fileInfo;
	string str;

	//check parameters
	if (k_mode!=1 && k_mode!=2) {cout<<"error! unknown kinship/relatedness input mode: "<<k_mode<<endl; error=true;}

	if (a_mode!=11 && a_mode!=12 && a_mode!=21 && a_mode!=31 && a_mode!=22 && a_mode!=43 && a_mode!=51 && a_mode!=52 && a_mode!=53 && a_mode!=54 && !saveGeno && !saveSS)
	{cout<<"error! unknown analysis mode: "<<a_mode<<". make sure -saveSS -saveGenoe -gk or -lm or -bvsrm or -predict is sepcified correctly."<<endl; error=true;}

	if (miss_level>1) {cout<<"error! missing level needs to be between 0 and 1. current value = "<<miss_level<<endl; error=true;}
	if (maf_level>0.5) {cout<<"error! maf level needs to be between 0 and 0.5. current value = "<<maf_level<<endl; error=true;}
	if (hwe_level>1) {cout<<"error! hwe level needs to be between 0 and 1. current value = "<<hwe_level<<endl; error=true;}

	if (h_max<h_min) {cout<<"error! maximum h value must be larger than the minimal value. current values = "<<h_max<<" and "<<h_min<<endl; error=true;}
	if (s_max<s_min) {cout<<"error! maximum s value must be larger than the minimal value. current values = "<<s_max<<" and "<<s_min<<endl; error=true;}
	if (rho_max<rho_min) {cout<<"error! maximum rho value must be larger than the minimal value. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
	if (logp_max<logp_min) {cout<<"error! maximum logp value must be larger than the minimal value. current values = "<<logp_max/log(10)<<" and "<<logp_min/log(10)<<endl; error=true;}

	if (h_max>1) {cout<<"error! h values must be bewtween 0 and 1. current values = "<<h_max<<" and "<<h_min<<endl; error=true;}
	if (rho_max>1) {cout<<"error! rho values must be between 0 and 1. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
	if (logp_max>0) {cout<<"error! maximum logp value must be smaller than 0. current values = "<<logp_max/log(10)<<" and "<<logp_min/log(10)<<endl; error=true;}

	if (h_scale>1.0) {cout<<"error! hscale value must be between 0 and 1. current value = "<<h_scale<<endl; error=true;}
	if (rho_scale>1.0) {cout<<"error! rscale value must be between 0 and 1. current value = "<<rho_scale<<endl; error=true;}
	if (logp_scale>1.0) {cout<<"error! pscale value must be between 0 and 1. current value = "<<logp_scale<<endl; error=true;}

	if (rho_max==1 && rho_min==1 && a_mode==12) {cout<<"error! ridge regression does not support a rho parameter. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
	// if (final_EM){cout<< "Final EM, calculate VB R2" << endl;}
	// if (!final_EM){cout<< "Not final EM, check " << endl;}

	//check if files are compatible with each other, and if files exist
	if (!file_bfile.empty()) {
		str=file_bfile+".bim";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .bim file: "<<str<<endl; error=true;}
		str=file_bfile+".bed";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .bed file: "<<str<<endl; error=true;}
		str=file_bfile+".fam";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .fam file: "<<str<<endl; error=true;}
	}


	str=file_geno;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open mean genotype file: "<<str<<endl; error=true;}

	size_t flag=0;
	if (!file_bfile.empty()) {flag++;}
	if (!file_geno.empty()) {flag++;}

	if (file_pheno.empty() && (a_mode==43) ) {
		cout<<"error! phenotype file is required."<<endl; error=true;
	}

	str=file_snps;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open snps file: "<<str<<endl; error=true;}


	str=file_anno;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open annotation file: "<<str<<endl; error=true;}

	str=file_kin;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open relatedness matrix file: "<<str<<endl; error=true;}

	//check if files are compatible with analysis mode
	if ((a_mode==43) && file_kin.empty())  {cout<<"error! missing relatedness file. -predict option requires -k option to provide a relatedness file."<<endl;  error=true;}

	if( inputSS && file_cov.empty() && file_score.empty() ) {
		cout << "Error! missing summary statistics file." << endl;
		error = true;
	}

	if( inputSS && (ni_test == 0 || rv == 0.0) ) {
		cout << "Error! -inputSS is specified, need input for -n [sample size] and -rv [phenotype variance] \n";
		error = true ;
	}

	return;
}




void PARAM::CheckData (void) {

	//calculate ni_total and ni_test, and set indicator_idv to 0 whenever indicator_cvt=0
	//and calculate np_obs and np_miss
	if(!inputSS){
		ni_total=(indicator_idv).size();

		ni_test=0;
		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			if (indicator_idv[i]==0) {continue;}
			ni_test++;
		}

		np_obs=0; np_miss=0;
		for (size_t j=0; j<indicator_pheno.size(); j++) {
				if (indicator_pheno[j]==0) {
					np_miss++;
				} else {
					np_obs++;
				}
		}

		if (ni_test==0) {
			error=true;
			cout<<"error! number of analyzed individuals equals 0. "<<endl;
			return;
		}
		if (a_mode==43) {
			cout<<"\n## number of observed data = "<<np_obs<<endl;
			cout<<"## number of missing data = "<<np_miss<<endl;
		}
		CreateSnpPosVec(snp_pos, snpInfo, indicator_snp);
    	// order snp_pos by chr/bp
    	stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp);
    	cout << "snp_pos size : " << snp_pos.size() << endl;
	}else{
		ni_total = ni_test;
		pheno_var = rv;
		// Create a vector of "SNPPOS" structs snp_pos (snpInfo will be cleared)
    	CreateSnpPosVec(snp_pos, snpInfo, indicator_snp);
    	// order snp_pos by chr/bp
    	stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp);
    	// Delete variants without LD information
    	UpdateScore();
    	// cout << "snp_pos size : " << snp_pos.size() << endl;
	}

	//output some information
	cout<<"\n## number of total individuals = "<<ni_total<<endl;
	cout<<"## number of individuals with full phenotypes = "<<ni_test<<endl;
	cout<<"## number of total SNPs = "<<ns_total<<endl;
	cout<<"## number of analyzed SNPs = "<<ns_test<<endl;

	//set parameters for BSLMM
	//and check for predict
	if (a_mode==11 || a_mode==12 || a_mode==13 || a_mode==31) {
		if (a_mode==11) {n_mh=1;}
		if (logp_min==0) {logp_min=-1.0*log((double)ns_test);}

		if (h_scale==-1) {h_scale=min(1.0, 10.0/sqrt((double)ni_test) );}
		if (rho_scale==-1) {rho_scale=min(1.0, 10.0/sqrt((double)ni_test) );}
		if (logp_scale==-1) {logp_scale=min(1.0, 5.0/sqrt((double)ni_test) );}
        //cout << "h_scale = " << h_scale << "; rho_scale = " << rho_scale<< "; logtheta_scale = " << logp_scale << endl;
        if (vscale <= 0.0) { vscale = min(0.5, 10.0/sqrt((double)ni_test));}

		if (h_min==-1) {h_min=0.00000001;}
		if (h_max==-1) {h_max=1.0;}

		if (s_max>ns_test) {s_max=ns_test; cout<<"s_max is re-set to the number of analyzed SNPs."<<endl;}
		if (s_max<s_min) {cout<<"error! maximum s value must be larger than the minimal value. current values = "<<s_max<<" and "<<s_min<<endl; error=true;}
	}

	return;
}


void PARAM::PrintSummary ()
{
	cout<<"pve estimate ="<<endl;
	cout<<"se(pve) ="<<endl;

	return;
}


void PARAM::ReadGenotypes (uchar **X, gsl_matrix *K) {

 	cout << "\nStarting reading genotype files for the second time ...\n";
    string file_str;
    UnCompBufferSize = (ni_test) * sizeof(uchar);
   // cout << "UnCompBufferSize = " << UnCompBufferSize << endl;

    if(calc_K) cout << "Kinship matrix is calculated here.\n";

	if (!file_bfile.empty()) {
		file_str=file_bfile+".bed";
		if (ReadFile_bed (file_str, indicator_idv, indicator_snp, X, K, calc_K, ni_test, ns_test, ni_total, ns_total, SNPmean, CompBuffSizeVec, Compress_Flag)==false) {error=true;}
        //revised
	}

    else if(!file_vcf.empty()){
        if ( ReadFile_vcf (file_vcf, indicator_idv, indicator_snp, X, ni_test, ns_test, K, calc_K, GTfield, SNPmean, CompBuffSizeVec, SampleVcfPos, PhenoID2Pos, VcfSampleID, Compress_Flag)==false )
        {error=true;} // revised
    }

    else if(!file_geno.empty()){
        if (ReadFile_geno (file_geno, indicator_idv, indicator_snp, X, K, calc_K, ni_test, SNPmean, CompBuffSizeVec, SampleVcfPos, PhenoID2Pos, VcfSampleID, Compress_Flag)==false) {error=true;} //to be revised
    }else{
    	cerr << "one of the genotype files has to be specified." << endl;
    	exit(-1);
    }

    return;

}

// Read summary statistics and load into cPar
void PARAM::ReadSS (){
	if( (! file_score.empty()) && (!file_cov.empty()) ){
    		cout << "\nLoad summary statistics ...\n";
    		if(ReadFile_corr(file_cov, ns_test, snpInfo, LD_ref, mapLDKey2Pos) == false)
    			{ error = true; }

    		if(ReadFile_score(file_score, snpInfo, mapScoreKey2Pos, mapLDKey2Pos, pval_vec, pos_ChisqTest, U_STAT, SQRT_V_STAT, xtx_vec, snp_var_vec, ns_test, ns_total, mbeta, mbeta_SE, indicator_snp, ni_test, maf_level, hwe_level, pheno_var, LD_ref, use_xtx_LD) == false)
    			{ error = true; }

    		// read functional/annotation fule
			if ( (!file_anno.empty()) && (!file_func_code.empty()) ) {
		    	cout << "\nStart loading annotation files ... \n";
		    	//cout << file_anno << " \nwith code file " << file_func_code << "\n";
		        if (ReadFile_anno (file_anno, file_func_code, mapScoreKey2Pos, mapFunc2Code, snpInfo, n_type, mFunc)==false)
		        	{error=true;}
		    }
        // remove new function for testing 11-17-18
		    else {
               cout << "\nAnnotation assignment by position!\n";
               if (setANNOCode (start_pos, end_pos, target_chr, window_size, indicator_snp, snpInfo, n_type, mFunc)==false) {error=true;}
            }
        //else {
            //if (Empty_anno (indicator_snp, snpInfo, n_type, mFunc)==false)
            //{error=true;}
        //}
    }else{
    	cerr << "Need to specify summary score.txt and cov.txt files!";
    	error = true;
    }
	return;
}


void PARAM::WriteGenotypes(uchar **X){

	string file_str;
    file_str="./output/"+file_out;
    file_str+=".geno";

    //cout << "create UcharTable ...\n";
    CreateUcharTable(UcharTable);

    ofstream outfile (file_str.c_str(), ofstream::out);
    if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

    //write header with VcfSampleID_test
    cout << "Write genotype dosage file for analyzed samples:" << ni_test << endl;
    //cout << "VcfSampleID_test length = " << VcfSampleID_test.size() << endl;

    outfile<<"#CHROM"<<"\t"<<"POS"<<"\t" <<"ID"<<"\t" << "REF"<< "\t" << "ALT"  << "\t";

    for (size_t i=0; i<ni_test; i++) {
    	//if(i < 10){ cout << VcfSampleID_test[i] <<endl; }
        if (i ==(ni_test-1)) {
            outfile << VcfSampleID_test[i] << endl;
        }
        else outfile << VcfSampleID_test[i] << "\t";
    }

    size_t pos=0;
    double geno_j;
    uchar c;

    //cout << "write variant information."<<endl;
    for (size_t i=0; i<ns_total; ++i) {

    	if(!indicator_snp[i]){continue;}

    	// save the data
        outfile<< snpInfo[i].chr<<"\t" <<snpInfo[i].base_position <<"\t"  << snpInfo[i].rs_number << "\t" << snpInfo[i].a_major << "\t" << snpInfo[i].a_minor << "\t";

        for (size_t j=0; j < ni_test; j++) {
        	c = X[pos][j];
            geno_j = UcharTable[(int)c].second;
            if(geno_j < 0.0 || geno_j > 2){
            	cout << "ERROR: genotype = " << geno_j <<" for " << snpInfo[i].rs_number <<":"<< snpInfo[i].chr<<":" <<snpInfo[i].base_position << ":" << snpInfo[i].a_major << ":" << snpInfo[i].a_minor << "; sample " << VcfSampleID_test[j] << endl;
                exit(-1);
            }else{
            		if (geno_j == 0.0) geno_j = 0;
            		else if (geno_j == 2.0) geno_j = 2;
            		else if (geno_j == 1.0) geno_j = 1;
		            if (j == (ni_test-1))
		                outfile << fixed << setprecision(2)  << geno_j << endl;
		            else
		                outfile << fixed << setprecision(2) << geno_j << "\t";
		            }
        }
        pos++;
    }

    outfile.clear();
    outfile.close();

}


// Calculate kinshiip matrix
void PARAM::CalcKin (gsl_matrix *matrix_kin)  {
	string file_str;

	gsl_matrix_set_zero (matrix_kin);

	if ( !file_bfile.empty() ) {
		file_str=file_bfile+".bed";
		if (PlinkKin (file_str, indicator_idv, indicator_snp, a_mode-20, d_pace, matrix_kin)==false) {error=true;}
	}
	else if( !file_geno.empty() ) {
		file_str=file_geno;
		if (GenoKin (file_str, indicator_idv, indicator_snp, a_mode-20, d_pace, matrix_kin, SampleVcfPos, PhenoID2Pos, VcfSampleID)==false) {error=true;}
	}
	else if( !file_vcf.empty() ) {
		file_str=file_vcf;
		if (VCFKin (file_str, indicator_idv, indicator_snp, a_mode-20, d_pace, matrix_kin, GTfield, SampleVcfPos, PhenoID2Pos, VcfSampleID)==false) {error=true;}
	}
	else{
		cerr << "Need input genotype file: plink, bimbam, or VCF!" <<endl;
	}

	return;
}





/*
void PARAM::CheckCvt ()
{
	if (indicator_cvt.size()==0) {return;}

	size_t ci_test=0;

	gsl_matrix *W=gsl_matrix_alloc (ni_test, n_cvt);

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0 || indicator_cvt[i]==0) {continue;}
		for (size_t j=0; j<n_cvt; ++j) {
			gsl_matrix_set (W, ci_test, j, (cvt)[i][j]);
		}
		ci_test++;
	}

	size_t flag_ipt=0;
	double v_min, v_max;
	set<size_t> set_remove;

	//check if any columns is an intercept
	for (size_t i=0; i<W->size2; i++) {
		gsl_vector_view w_col=gsl_matrix_column (W, i);
		gsl_vector_minmax (&w_col.vector, &v_min, &v_max);
		if (v_min==v_max) {flag_ipt=1; set_remove.insert (i);}
	}

	//add an intecept term if needed
	if (n_cvt==set_remove.size()) {
		indicator_cvt.clear();
		n_cvt=1;
	} else if (flag_ipt==0) {
		cout<<"no intecept term is found in the cvt file. a column of 1s is added."<<endl;
		for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
			if (indicator_idv[i]==0 || indicator_cvt[i]==0) {continue;}
			cvt[i].push_back(1.0);
		}

		n_cvt++;
	} else {}

	gsl_matrix_free(W);

	return;
}
*/

//reorder phenotypes
void PARAM::ReorderPheno(gsl_vector *y)
{
	cout << "Reorder phenotype with respect to sample orders in the vcf/geno files... "<< endl;

    double pheno_i;
    string id;
    size_t c_ind=0, pheno_idx;
    gsl_vector *ytemp=gsl_vector_alloc (ni_test);
    VcfSampleID_test.clear(); // set VCFSampleID_test

    cout << "Total number of samples in the genotype file: " << VcfSampleID.size() << endl;
	cout << "Total number of samples in the phenotype file: " << PhenoID2Pos.size() << endl;
	cout << "Total number of samples in with both geno and pheno data: " << PhenoID2Ind.size() << endl;

	// indicator_idv is of the same order as in the phenotype file
	GenoSampleID2Ind.clear();

    for (size_t i=0; i < VcfSampleID.size(); i++) {
        id = VcfSampleID[i];

        if( (PhenoID2Ind.count(id) > 0) & (GenoSampleID2Ind.count(id) == 0) ){
        	//if (i < 10) cout << id << ", count ="<< PhenoID2Ind.count(id) << ";    ";
        	pheno_idx = PhenoID2Ind[id];

        	pheno_i = gsl_vector_get(y, pheno_idx);
            gsl_vector_set(ytemp, c_ind, pheno_i);
            VcfSampleID_test.push_back(id);
            GenoSampleID2Ind[id] = i;
            c_ind++;

        }
    }
   // cout << "Finished reordering phenotypes. \nAnalyzed sample size is " << c_ind << endl;
   // cout << "GenoSampleID2Ind.size() = " << GenoSampleID2Ind.size() << endl;
    ni_test = c_ind;
    gsl_vector_memcpy(y, ytemp);
    gsl_vector_free(ytemp);
}


// Post-process phentoypes, obtain ni_total, ni_test, PhenoID2Ind
void PARAM::ProcessPheno()
{
	//obtain ni_test
	ni_test=0;

	//obtain PhenoID2Ind map (will need to adjust for samples with geno)
    PhenoID2Ind.clear(); // map sample ID to index in analyzed pheno
    PhenoID2Pos.clear(); // map sample ID to index in raw/all pheno

	for (size_t i=0; i<ni_total; ++i) {
		PhenoID2Pos[InputSampleID[i]] = i;

        if (indicator_idv[i]){
            PhenoID2Ind[InputSampleID[i]]= ni_test;
            ni_test++;
        }
      //  else PhenoID2Ind[InputSampleID[i]] = ULONG_MAX;
    }
    //cout << "Create PhenoID2Ind map; total number of analyzed individual ni_test = " << ni_test << "\n";
    //cout << "PhenoID2Ind map length = " << PhenoID2Ind.size() << "\n";

	if (ni_test==0) {
		error=true;
		cout<<"error! number of analyzed individuals equals 0. "<<endl;
		return;
	}

	return;
}

// Update phentoype sample info: ni_test, PhenoID2Ind
void PARAM::UpdatePheno()
{
	//cout << "Update ni_test, PhenoID2Ind for samples that have geno data based on GenoSampleID2Ind ...\n";
	//cout << "GenoSampleID2Ind map length = " << GenoSampleID2Ind.size() << "\n";

	//obtain ni_test
	ni_test=0;

	//obtain PhenoID2Ind map (will need to adjust for samples with geno)
    PhenoID2Ind.clear();

	for (size_t i=0; i<ni_total; ++i) {
        if (indicator_idv[i]){
        	if( GenoSampleID2Ind.count(InputSampleID[i]) > 0 ){
        		PhenoID2Ind[InputSampleID[i]]= ni_test;
            	ni_test++;
        	}else{
        		indicator_idv[i] = false;
        	}
        }
      //  else PhenoID2Ind[InputSampleID[i]] = ULONG_MAX;
    }
    // cout << "Update PhenoID2Ind map; \n";
    cout << "Total number of analyzed individual with geno: ni_test = " << ni_test << "\n";
    cout << "PhenoID2Ind map length = " << PhenoID2Ind.size() << "\n";

	if (ni_test==0) {
		error=true;
		cout<<"error! number of analyzed individuals equals 0. "<<endl;
		return;
	}

	return;
}



//else, indicator_idv to load phenotype
void PARAM::CopyPheno (gsl_vector *y)
{
	size_t ci_test=0;
	pheno_mean = 0.0;

	for (size_t i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]) {
			gsl_vector_set (y, ci_test, pheno[i]);
			pheno_mean += pheno[i];
			ci_test++;
		}
	}
	pheno_mean /= (double) ci_test;
	gsl_vector_add_constant(y, -1.0 * pheno_mean);

	return;
}


void CreateSnpPosVec(vector<SNPPOS> &snp_pos, vector<SNPINFO> &snpInfo, const vector<bool> &indicator_snp)
{
    size_t pos;
    string rs;
    string chr;
    long int bp;
    size_t tt=0;
    vector<bool> indicator_func;
    //vector<double> weight;
    //double weight_i;
    double maf;
    string a_minor;
    string a_major;
    string key;

    snp_pos.clear();

    for (size_t i=0; i < indicator_snp.size(); ++i){

        if(!indicator_snp[i]) {continue;}

        pos=tt;
       // if(tt == pos_loglr[tt].first ) pos = tt;
        //else cout << "error assigning position to snp_pos vector"<< endl;

        rs = snpInfo[i].rs_number;
        chr = snpInfo[i].chr;
        bp = snpInfo[i].base_position;
        maf = snpInfo[i].maf;
        a_minor = snpInfo[i].a_minor;
        a_major = snpInfo[i].a_major;
        key = snpInfo[i].key;

        indicator_func = snpInfo[i].indicator_func;

        SNPPOS snp_temp={pos, rs, chr, bp, a_minor, a_major, maf, indicator_func, key};
        snp_pos.push_back(snp_temp);

        tt++;
    }
    snpInfo.clear();

}


// Update variants in the summary stat vectors
void PARAM::UpdateScore(){
	cout << "\nns_test from summary stat is " << ns_test << endl;
    double beta2_i, beta_se2_i, u_i, v_i, ni_effect_i, xtx_i;
    double yty = pheno_var * ((double)ni_test - 1.0);
    // cout << "yty in UpdateScore is " << yty << endl;

    ni_effect_vec.clear();
    if(scaleN){
    	cout << "Scale summary statistics with effective sample size!\n";
    	for(size_t i=0; i<ns_test; i++){
			// Scale xtx_i value, update score statistics
			beta2_i = mbeta[i] * mbeta[i];
			beta_se2_i = mbeta_SE[i] * mbeta_SE[i];
			ni_effect_i = ni_test;
			xtx_i = yty / ((beta_se2_i * (double)(ni_test - 1) ) + beta2_i) ;
			// ni_effect_i = (yty / (xtx_vec[i] * beta_se2_i)) - (beta2_i / beta_se2_i) + 1;
			//if(ni_effect_i < 2) ni_effect_i = 2;
			//xtx_i = (double)(ni_effect_i) * snp_var_vec[i];
			u_i = xtx_i * mbeta[i] ;
			v_i = sqrt(xtx_i * pheno_var);
			ni_effect_vec.push_back(ni_effect_i);
			U_STAT[i] = u_i;
			SQRT_V_STAT[i] = v_i;
			xtx_vec[i] = xtx_i;
		}
    }else{
    	ni_effect_vec.assign(ns_test, ni_test);
    }

	// Calculate trace of X
	trace_G=0;
	for(size_t i = 0; i < xtx_vec.size() ; i++){
		trace_G +=  xtx_vec[i];
	}
}

// Get LD matrix from Ref LDR2.txt for variants in score.txt file
void PARAM::Convert_LD(){

// convert cov matrix to LD r2 matrix, using previously saved xtx vector
    double r2_ij;
    string key_i, key_j;
    size_t idx_i, idx_j;
    bool swap_i, swap_j;

    // cout << "Convert_LD ns_test = " << ns_test << "; ni_test = " << ni_test << endl;
    LD.clear();

    // Set up LD matrix for all variants in score.txt
    for(size_t i=0; i<snp_pos.size(); i++){

        LD.push_back(vector<double>());
        LD[i].push_back(snp_var_vec[i]) ;
        key_i = snp_pos[i].key;
        swap_i = false;

        if( mapLDKey2Pos.count(key_i) == 0 ){
            SwapKey(key_i);
            swap_i = true;
            if( mapLDKey2Pos.count(key_i) == 0 ){
                cout << key_i << " dose not appear in cov.txt file\n";
                exit(-1);
            }
            else{
                idx_i = mapLDKey2Pos[key_i];
                for(size_t j = (i+1); j < snp_pos.size(); j++)
                {
                    key_j = snp_pos[j].key;
                    swap_j = false;
                    if( mapLDKey2Pos.count(key_j) == 0 ){
                        SwapKey(key_j);
                        if( mapLDKey2Pos.count(key_j) == 0 ){
                            cout << key_j << " dose not appear in cov.txt file\n";
                            exit(-1);
                        }else{
                            swap_j = true;
                            idx_j = mapLDKey2Pos[key_j];
                            if( abs( (int)idx_j - (int)idx_i) >= LD_ref[ min(idx_i, idx_j) ].size() )
                        		continue;
                            r2_ij = getR2_ij(LD_ref, idx_i, idx_j, swap_i, swap_j);
                            if(refLD & (r2_ij < r2_level) )
                        		r2_ij = 0.0;
                            LD[i].push_back(r2_ij);
                        }
                    }else{
                        idx_j = mapLDKey2Pos[key_j];
                        if( abs( (int)idx_j - (int)idx_i) >= LD_ref[ min(idx_i, idx_j) ].size() )
                        		continue;
                        r2_ij = getR2_ij(LD_ref, idx_i, idx_j, swap_i, swap_j);
                        if(refLD & (r2_ij < r2_level))
                        	r2_ij = 0.0;
                        LD[i].push_back(r2_ij);
                    }
                }
            }
        }else{
            idx_i = mapLDKey2Pos[key_i];
            for(size_t j = (i+1); j < snp_pos.size(); j++)
            {
                key_j = snp_pos[j].key;
                swap_j = false;
                if( mapLDKey2Pos.count(key_j) == 0 ){
                    SwapKey(key_j);
                    swap_j = true;
                    if( mapLDKey2Pos.count(key_j) == 0 ){
                        cout << key_j << " dose not appear in cov.txt file\n";
                        exit(-1);
                    }else{
                        idx_j = mapLDKey2Pos[key_j];
                        if( abs((int)idx_j - (int)idx_i) >= LD_ref[ min(idx_i, idx_j) ].size() )
                        		continue;
                        r2_ij = getR2_ij(LD_ref, idx_i, idx_j, swap_i, swap_j);
                        if(refLD & (r2_ij < r2_level))
                        	r2_ij = 0.0;
                        LD[i].push_back(r2_ij);
                    }
                }else{
                    idx_j = mapLDKey2Pos[key_j];
                    if( abs((int)idx_j - (int)idx_i) >= LD_ref[ min(idx_i, idx_j) ].size() )
                        		continue;
                    r2_ij = getR2_ij(LD_ref, idx_i, idx_j, swap_i, swap_j);
                    if(refLD & (r2_ij < r2_level))
                        	r2_ij = 0.0;
                    LD[i].push_back(r2_ij);
                }
            }
        }
    }
    // clear LD_ref to save memory
    LD_ref.clear();
    return;
}



vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}


void SwapKey(string &key){
	vector<string> temp;
	temp = split(key, ":");
	if(temp.size() < 4) {
		// cout << "temp.size" <<  temp.size << endl;
		cout << "SNP annotation (chr:pos:ref:alt) is in wrong format: " << key << endl;
		exit(-1);
	}else{
		key = temp[0] + ":" + temp[1] + ":" + temp[3] + ":" + temp[2] ;
	}
	return;
}

// Obtain correlation between snp_i and snp_j
double getR2_ij(const vector< vector<double> > &LD, size_t idx_i, size_t idx_j, const bool &swap_i, const bool &swap_j)
{

    double r2_ij = 0.0;
    size_t idx_temp, idx_dist;

    // make pos_i the front position
    if(idx_i > idx_j){
        idx_temp = idx_i;
        idx_i = idx_j;
        idx_j = idx_temp;
    }

    idx_dist = idx_j - idx_i;
    if(idx_dist == 0){
        r2_ij = 1.0;
    }
    else if( idx_dist < LD[idx_i].size() )
    {
        r2_ij = LD[idx_i][idx_dist];
        if( (swap_i + swap_j) == 1 ){
            r2_ij = -r2_ij ;
        }
    }
    return r2_ij;
}
