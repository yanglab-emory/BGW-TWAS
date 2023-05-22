options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)

gene_name <- args[[1]]
wkdir <- args[[2]]

print( paste0("gene name: ", gene_name))
print( paste0("Working Directory: ", wkdir))

library(scales)
library(data.table)
setwd(paste0(wkdir))

# gene_name="CCDC106"
### Load pred_grex.txt and weight files
grex_pred <- read.table(paste0(gene_name, "_pred_grex.txt"), header=TRUE)
rownames(grex_pred) <- grex_pred$sampleID

ge_dt <- fread(paste0(gene_name, "_ge_train.txt" ), header=TRUE)
ge_dt <- as.data.frame(ge_dt)
# print(ge_dt[, 1:5])
ge_vec <- unlist(ge_dt[1, grex_pred$sampleID])

sum.fit <- summary(lm(ge_vec ~ grex_pred[, 2]))
pvalue <- pf(sum.fit$fstatistic[1], sum.fit$fstatistic[2], sum.fit$fstatistic[3], lower.tail=FALSE)


print(data.frame(Gene=gene_name, TrainR2=round(sum.fit$r.squared, 3), PVALUE=scientific(pvalue, digits=3)))

write.table(data.frame(Gene=gene_name, TrainR2=round(sum.fit$r.squared, 3), PVALUE=scientific(pvalue, digits=3)),
	file=paste0(gene_name, "_TrainR2.txt"),
	col.names=TRUE, row.names=FALSE, quote=FALSE)

