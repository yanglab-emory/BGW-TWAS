options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)

gene_name<-args[1]
wkdir<-args[2]

library(data.table)
setwd(paste0(wkdir))

# gene_name="ABCA7"
### Load grex.geno and weight files
geno_pred<- fread(paste0(gene_name, "_grex.geno"), header=TRUE)
geno_pred$ID <- paste(geno_pred$CHROM, geno_pred$POS, sep = ":")
setkey(geno_pred, "ID")

BGW_weights <- fread(paste0(gene_name, "_BGW_eQTL_weights.txt" ), header=TRUE)
BGW_weights <- BGW_weights[, c("CHROM", "POS", "ID", "REF", "ALT",  "Trans", "PCP", "beta")]
BGW_weights$ID <- paste(BGW_weights$CHROM, BGW_weights$POS, sep = ":")
BGW_weights$w = BGW_weights$PCP * BGW_weights$beta
setkey(BGW_weights, "ID")

### Calculate GReX, sum cis PCP, sum trans PCP
tmp_test = NULL; w_vec = NULL; cis_PCP = 0; trans_PCP = 0;
n_snp = nrow(BGW_weights)
if( n_snp > 0 ){

	for(i in 1:n_snp){
		temp_snp_id = BGW_weights[i, "ID"]
		ref = BGW_weights[i, "REF"]
		alt = BGW_weights[i, "ALT"]

		temp_geno = geno_pred[BGW_weights[i, "ID"], ]
		if(!is.na(temp_geno$POS) & temp_geno$REF == ref & temp_geno$ALT == alt){
			tmp_test = rbind(tmp_test, c(BGW_weights[i, c("CHROM", "POS", "ID", "REF", "ALT")], geno_pred[BGW_weights[i, "ID"], -c(1:5)]) )
			w_vec = c(w_vec, BGW_weights[i, ]$w)
			if(BGW_weights[i, ]$Trans){
				trans_PCP = trans_PCP + BGW_weights[i, ]$PCP
			}else{
				cis_PCP = cis_PCP + BGW_weights[i, ]$PCP
			}
		}
		else if(!is.na(temp_geno$POS) & temp_geno$REF == alt & temp_geno$ALT == ref){
			tmp_test = rbind(tmp_test, c(BGW_weights[i, c("CHROM", "POS", "ID", "REF", "ALT")], 2 - geno_pred[BGW_weights[i, "ID"], -c(1:5)]) )
			w_vec = c(w_vec, BGW_weights[i, ]$w)
			if(BGW_weights[i, ]$Trans){
				trans_PCP = trans_PCP + BGW_weights[i, ]$PCP
			}else{
				cis_PCP = cis_PCP + BGW_weights[i, ]$PCP
			}
		}
		else{
			print(paste("SNP", temp_snp_id, ref, alt, "do not have test genotype!" ))
		}
	}
}else{
	print(gene_name, "has no eQTL in the weight file!")
}


# names(gwas_res)<- c("CHROM", "POS", "ID", "REF", "ALT", "maf", "location", "pi", "beta")

X_geno <- matrix( as.numeric(t(tmp_test[, -(1:5)]) ), ncol = length(w_vec))
is.numeric(X_geno)
GREX_hat <- X_geno %*% (w_vec)
write.table(data.frame(sampleID=colnames(tmp_test)[-(1:5)], GReX=GREX_hat)
,
file=paste0(gene_name,"_pred_grex.txt"),
quote = FALSE, sep ="\t", row.names = FALSE)


total_PCP<-cis_PCP+trans_PCP
write.table(data.frame(Gene=gene_name, total_PCP=total_PCP, cis_PCP=cis_PCP, trans_PCP=trans_PCP),
	file=paste0(gene_name, "_sumPCP.txt"),
	col.names=TRUE, row.names=FALSE, quote=FALSE)

