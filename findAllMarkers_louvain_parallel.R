require(Seurat); require(dplyr)

load("/home/dcha/02.glioblastoma_scRNAseq/rdata/obj.fastmnn_pc20.RData")
DefaultAssay(obj.fastmnn_pc20)="RNA"
Idents(obj.fastmnn_pc20)=obj.fastmnn_pc20$seurat_clusters

args=commandArgs(trailingOnly = TRUE)

testUse=args[1]
clusterOfInterest=args[2]
outputFn=paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/fastMNN_markers/cluster_",clusterOfInterest, "_markers.",testUse, ".rds")
print(paste0(testUse," for detect markers of ", clusterOfInterest))
print(paste0("RDS: ",paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/fastMNN_markers/cluster_",clusterOfInterest, "_markers.",testUse, ".rds")))

others=Idents(obj.fastmnn_pc20) %>% levels(.) %>% setdiff(., clusterOfInterest)

clustermarkers= FindMarkers(object = obj.fastmnn_pc20, ident.1=clusterOfInterest, ident.2=others, 
                            test.use=testUse, 
                            only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(clustermarkers, file=outputFn)
