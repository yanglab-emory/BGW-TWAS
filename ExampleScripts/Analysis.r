library(tidyverse)
library(gridExtra)
source("~/GIT/BGW-TWAS/bin/R_funcs.r")

wkdir= "/mnt/YangFSS/data2/jyang/testBGW/TP63_wk/" # BGW-TWAS EM_MCMC directory
setwd(wkdir)

########
bgw_weights = fread("./TP63_BGW_eQTL_weights.txt", header = TRUE)
head(bgw_weights)

## Number of expected eQTL: Sum of PCP
sum(bgw_weights$PCP)

#### n_eQTL: Numbers of expected eQTL for cis- and trans- groups
#### beta_var: Effect size variance for cis- and trans- eQTL
bgw_weights %>% group_by(Trans) %>% summarise(n_eQTL = sum(PCP), beta_var = var(beta))

## Likely significant causal eQTL
bgw_weights[bgw_weights$PCP > 0.1068, ]

###### Plot eQTL results by BGW-TWAS
plot_pcp <- myManPlot_PCP(manPlot_dt = bgw_weights, chr_vec = sort(unique(bgw_weights$`#CHROM`)),
	title = "eQTL PCP by BGW-TWAS", chrGAP = 10)

plot_beta <- myManPlot_beta(manPlot_dt = bgw_weights, chr_vec = sort(unique(bgw_weights$`#CHROM`)),
	title = "eQTL effect size by BGW-TWAS", chrGAP = 10)

plot_weight <- myManPlot_weight(manPlot_dt = bgw_weights, chr_vec = sort(unique(bgw_weights$`#CHROM`)),
	title = "eQTL weight by BGW-TWAS", chrGAP = 10)

pdf("/mnt/YangFSS/data2/jyang/testBGW/TP63_wk/AnalysisResults/TP63_eQTL.pdf", width = 12, height = 10)
grid.arrange(plot_pcp , plot_beta , plot_weight, ncol = 1, nrow = 3)
dev.off()















