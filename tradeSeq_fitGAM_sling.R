#!/usr/bin/Rscript

require(magrittr); require(tradeSeq); require(slingshot); require(dplyr)

args=commandArgs(trailingOnly = TRUE)
fn=args[1]
path="/home/dcha/02.glioblastoma_scRNAseq/rdata/slingshot_data/"
sce=readRDS(paste0(path,fn))
fitGAM_sce = fitGAM(sce, nknots=5)

fitGAM_sce_fn = fn %>% gsub(".rds", ".fitGAM.rds", .)
saveRDS(fitGAM_sce, paste0(path, fitGAM_sce_fn))
