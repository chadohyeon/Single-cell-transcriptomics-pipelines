#!/usr/bin/Rscript
require(cFIT); require(dplyr); require(Seurat)

args=commandArgs(trailingOnly = TRUE)
toCFIT=args[1]
nCores=args[2]

cFitOut=readRDS(paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cfit_",toCFIT, ".rds"))
int.out=cFitOut$int.out
exprs.list=cFitOut$exprs.list
data.list=cFitOut$data.list

obj=readRDS("/home/dcha/02.glioblastoma_scRNAseq/rdata/merged_hart.rds")

mouseTerms=obj@meta.data %>% filter(species=="mouse") %>% .[[toCFIT]] %>% unique(.)
humanTerms=obj@meta.data %>% filter(species=="human") %>% .[[toCFIT]] %>% unique(.)

ref.data=which(names(exprs.list) %in% mouseTerms)
query.data=which(names(exprs.list) %in% humanTerms)

tf.out = CFITTransfer(Xtarget=exprs.list[[query.data]], Wref=int.out$W, max.niter = 100, seed=0, verbose=F, n.cores = nCores)

est.labels = asign_labels(exprs.source=do.call(rbind, int.out$H.list), 
                          exprs.target=tf.out$H, 
                          labels.source=do.call(c, data.list$labels.list[ref.data]))

cFitOut[["tf.out"]]=tf.out
cFitOut[["est.labels"]]=est.labels

saveRDS(cFitOut, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cfit_",toCFIT, ".rds"))