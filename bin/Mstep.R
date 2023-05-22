Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=F)
# source("/home/jyang/GIT/bfGWAS_SS/bin/R_funcs.r")

######## Need to pass args(hypfile, paramfile, k, hypcurrent_file) from bash
args <- commandArgs(TRUE)
# print(args)

## JML: add option to flag VB inference, and modify the data input by hypfile and EM result
## if using VB for obtaining transcriptome weights, we don't need the likelihood or the GVEC. we get different output from E step.

hypfile=args[[1]]
k=as.numeric(args[[2]])
pp = as.numeric(args[[3]])
abgamma = as.numeric(args[[4]])
EM_result_file = args[[5]]
hypcurrent_file=args[[6]]
#if (length(args) == 7) 
#VB_flag = 7
## need to ensure that the length of args is 6, not 7, if the 7th flag = "". 

## alternatively, if args[[7]] == "VB", do, else do something different. and update perl script to have a 7th flag that is empty or MCMC.  

if (args[[7]]=="VB"){ ## do VB version of Mstep
  
  # abgamma=0.1
  a_gamma <- b_gamma <- abgamma;
  print(paste("a_gamma=b_gamma = ", abgamma))
  print("updating with VB estimates")
  ###### Define functions to be used
  
  ## JML Change 3/18/19: in the VB function, we calculated the product phi_i*m_i^2 directly, 
  # rather than only m_i^2. do not need to weigh again by m in the sum. 
  CI_fish_pi <- function(m, p, a, b){
    pi_hat = (sum(m) + a - 1.0) / (p + a + b - 2.0)
    if(pi_hat <= 0 || pi_hat > 1){
      pi_hat = a/(a+b)
      se_pi = 0
    }else{
      if(pi_hat * (1-pi_hat) / (p + a + b - 2.0) > 0){
        se_pi = sqrt(pi_hat * (1-pi_hat) / (p + a + b - 2.0))
      }else{se_pi=0}
    }
    return(c(pi_hat, se_pi))
  }
  # in the VB function, we obtain directly sum over i phi*m_i*m_i, so we only need to sum the sigma2 term.
  Est_sigma2 <- function(sigma2, m, tau, a, b){
    sigma2_hat = (sum(sigma2) * tau + 2 * b) / (sum(m) + 2 * (a + 1))
    return(sigma2_hat)
  }
  
  CI_fish_sigma2 <- function(sigma2, m, tau, a, b){
    sigma2_hat = Est_sigma2(sigma2, m, tau, a, b)
    if( (sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat) < 0){
      se_sigma2=0
    }else{
      se_sigma2 = sigma2_hat * sqrt(1/(sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat))
    }
    return(c(sigma2_hat, se_sigma2))
  }
  
  ## log prior functions
  logprior_sigma <- function(a, b, x){ return(-(1+a) * log(x) - b/x) }
  
  logprior_pi <- function(a, b, x){ return((a-1) * log(x) + (b-1) * log(1 - x)) }
  
  
  ptm <- proc.time()
  ########### Load hypfile ....
  
  # hypfile="/net/fantasia/home/yjingj/GIT/bfGWAS/1KG_example/Test_Wkdir/hypval.current 
  
  hypdata = read.table(hypfile, sep="\t", header=FALSE)
  n_type = (dim(hypdata)[2] - 4)/4
  print(paste(" Total Annotation categories : ", n_type))
  
  temp_col_names <- c("block", "elbo_min", "GV", "rv")
  for(i in 1:n_type){
    temp_col_names <- c(temp_col_names, 
                        paste(c("n", "G", "m", "sigma2"), (i-1), sep = "_"))
  }
  
  colnames(hypdata) <-  temp_col_names
  print(head(hypdata)) 
  
  ########### Update hyper parameter values
  rv = mean(hypdata[, "rv"])
  tau = 1.0 / rv
  pve = sum(hypdata[, "GV"])
  
  prehyp <- read.table(hypcurrent_file, header=TRUE)
  print("hyper parameter values before VB: ")
  print(prehyp)
  
  ######### Set hierarchical parameter values
  n_vec = rep(0, n_type)
  for(i in 1:n_type){
    n_vec[i] <- sum(hypdata[, paste("n", (i-1), sep="_")])
  }
  print("checking n_snp counts: ")
  print(n_vec)
  #### updating hyper pi and sigma2 values for each group
  hypcurrent <- NULL
  hypmat <- NULL
  
  for(i in 1:n_type){
    # print(i)
    if(n_vec[i] > 0){
      a_beta = 2 * n_vec[i] * pp; b_beta = 2 * n_vec[i] - a_beta;
    }else{a_beta=1; b_beta = 1e6 - 1;}
    
    m_temp = hypdata[, paste("m", (i-1), sep="_")]
    
    pi_temp = CI_fish_pi(m_temp, n_vec[i], a_beta, b_beta)
    
    sigma2_temp = CI_fish_sigma2(hypdata[, paste("sigma2", (i-1), sep="_")], m_temp, tau, a_gamma, b_gamma)
    
    hypcurrent <- c(hypcurrent, pi_temp, sigma2_temp)
    # print(cbind(pi_temp, sigma2_temp))
    hypmat <- rbind(hypmat, c(pi_temp[1], sigma2_temp[1]))
  }
  
  ########## Write out updated hyper parameter values
  colnames(hypmat) <- c("pi", "sigma2")
  print("hyper parameter values updates after VB: ")
  print(hypmat)
  write.table(format(hypmat, scientific=TRUE), 
              file=hypcurrent_file, 
              quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
  
  #### Summarize log-likelihood
  loglike_total = sum(hypdata$loglike)
  
  for(i in 1:n_type){
    if(sum(prehyp[i, ]>0)==2){
      
      if(n_vec[i] > 0){
        a_beta = 2 * n_vec[i] * pp; b_beta = 2 * n_vec[i] - a_beta;
      }else{a_beta=1; b_beta = 1e6 - 1;}
      
      loglike_total = loglike_total + 
        logprior_pi(a_beta, b_beta, prehyp[i, 1]) +
        logprior_sigma(a_gamma, b_gamma, prehyp[i, 2]) 
    }else{
      print("pre-hyper-parameter <= 0... ")
    }
  }
  
  ########## Write out updated hyper parameter values and se to EM_result_file
  # EM_result_file="/net/fantasia/home/yjingj/GIT/bfGWAS/1KG_example/Test_Wkdir/Eoutput/EM_result.txt"
  hypcurrent = c(pve, loglike_total, hypcurrent)
  hypcurrent <- format(hypcurrent, scientific = TRUE)
  print("write to hypcurrent file with hyper parameter values after VB: ")
  print(c(k, hypcurrent))
  write.table(matrix(c(k, hypcurrent), nrow=1), file = EM_result_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
  
  print("EM step time cost (in minutes) : ")
  print((proc.time() - ptm)/60)

  ## end VB Mstep
  
} else { # abgamma=0.1
  a_gamma <- b_gamma <- abgamma;
  print(paste("a_gamma=b_gamma = ", abgamma))
  
  ###### Define functions to be used
  CI_fish_pi <- function(m, p, a, b){
    pi_hat = (sum(m) + a - 1.0) / (p + a + b - 2.0)
    if(pi_hat <= 0 || pi_hat > 1){
      pi_hat = a/(a+b)
      se_pi = 0
    }else{
      if(pi_hat * (1-pi_hat) / (p + a + b - 2.0) > 0){
        se_pi = sqrt(pi_hat * (1-pi_hat) / (p + a + b - 2.0))
      }else{se_pi=0}
    }
    return(c(pi_hat, se_pi))
  }
  
  Est_sigma2 <- function(sigma2, m, tau, a, b){
    sigma2_hat = (sum(sigma2 * m) * tau + 2 * b) / (sum(m) + 2 * (a + 1))
    return(sigma2_hat)
  }
  
  CI_fish_sigma2 <- function(sigma2, m, tau, a, b){
    sigma2_hat = Est_sigma2(sigma2, m, tau, a, b)
    if( (sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat) < 0){
      se_sigma2=0
    }else{
      se_sigma2 = sigma2_hat * sqrt(1/(sum(m) * tau - sum(m)/2 - (a+1) + 2*b/sigma2_hat))
    }
    return(c(sigma2_hat, se_sigma2))
  }
  
  ## log prior functions
  logprior_sigma <- function(a, b, x){ return(-(1+a) * log(x) - b/x) }
  
  logprior_pi <- function(a, b, x){ return((a-1) * log(x) + (b-1) * log(1 - x)) }
  
  
  ptm <- proc.time()
  ########### Load hypfile ....
  
  # hypfile="/net/fantasia/home/yjingj/GIT/bfGWAS/1KG_example/Test_Wkdir/hypval.current 
  
  hypdata = read.table(hypfile, sep="\t", header=FALSE)
  n_type = (dim(hypdata)[2] - 4)/4
  print(paste(" Total Annotation categories : ", n_type))
  
  temp_col_names <- c("block", "loglike", "GV", "rv")
  for(i in 1:n_type){
    temp_col_names <- c(temp_col_names, 
                        paste(c("n", "G", "m", "sigma2"), (i-1), sep = "_"))
  }
  
  colnames(hypdata) <-  temp_col_names
  
  ########### Update hyper parameter values
  rv = mean(hypdata[, "rv"])
  tau = 1.0 / rv
  pve = sum(hypdata[, "GV"])
  
  prehyp <- read.table(hypcurrent_file, header=TRUE)
  print("hyper parameter values before MCMC: ")
  print(prehyp)
  
  ######### Set hierarchical parameter values
  n_vec = rep(0, n_type)
  for(i in 1:n_type){
    n_vec[i] <- sum(hypdata[, paste("n", (i-1), sep="_")])
  }
  
  #### updating hyper pi and sigma2 values for each group
  hypcurrent <- NULL
  hypmat <- NULL
  
  for(i in 1:n_type){
    # print(i)
    if(n_vec[i] > 0){
      a_beta = 2 * n_vec[i] * pp; b_beta = 2 * n_vec[i] - a_beta;
    }else{a_beta=1; b_beta = 1e6 - 1;}
    
    m_temp = hypdata[, paste("m", (i-1), sep="_")]
    
    pi_temp = CI_fish_pi(m_temp, n_vec[i], a_beta, b_beta)
    
    sigma2_temp = CI_fish_sigma2(hypdata[, paste("sigma2", (i-1), sep="_")], m_temp, tau, a_gamma, b_gamma)
    
    hypcurrent <- c(hypcurrent, pi_temp, sigma2_temp)
    # print(cbind(pi_temp, sigma2_temp))
    hypmat <- rbind(hypmat, c(pi_temp[1], sigma2_temp[1]))
  }
  
  ########## Write out updated hyper parameter values
  colnames(hypmat) <- c("pi", "sigma2")
  print("hyper parameter values updates after MCMC: ")
  print(hypmat)
  write.table(format(hypmat, scientific=TRUE), 
              file=hypcurrent_file, 
              quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
  
  #### Summarize log-likelihood
  loglike_total = sum(hypdata$loglike)
  
  for(i in 1:n_type){
    if(sum(prehyp[i, ]>0)==2){
      
      if(n_vec[i] > 0){
        a_beta = 2 * n_vec[i] * pp; b_beta = 2 * n_vec[i] - a_beta;
      }else{a_beta=1; b_beta = 1e6 - 1;}
      
      loglike_total = loglike_total + 
        logprior_pi(a_beta, b_beta, prehyp[i, 1]) +
        logprior_sigma(a_gamma, b_gamma, prehyp[i, 2]) 
    }else{
      print("pre-hyper-parameter <= 0... ")
    }
  }
  
  ########## Write out updated hyper parameter values and se to EM_result_file
  # EM_result_file="/net/fantasia/home/yjingj/GIT/bfGWAS/1KG_example/Test_Wkdir/Eoutput/EM_result.txt"
  hypcurrent = c(pve, loglike_total, hypcurrent)
  hypcurrent <- format(hypcurrent, scientific = TRUE)
  print("write to hypcurrent file with hyper parameter values after MCMC: ")
  print(c(k, hypcurrent))
  write.table(matrix(c(k, hypcurrent), nrow=1), file = EM_result_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
  
  print("EM step time cost (in minutes) : ")
  print((proc.time() - ptm)/60)
  

}







