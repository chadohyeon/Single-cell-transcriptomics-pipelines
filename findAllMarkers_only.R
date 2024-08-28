#!/usr/bin/Rscript
require(Seurat); require(dplyr)

args=commandArgs(trailingOnly = TRUE)
rds=args[1]
obj=readRDS(paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/human_corrected/", rds))
Idents(obj)=obj$seurat_clusters

outputFn=paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/human_corrected/wilcoxRes.", rds)

clustermarkers=FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(clustermarkers, file=outputFn)
