dyn.load('/home/public/tools/miniconda3/lib/libhdf5.so.103')
require(rliger)

require(Seurat); require(SeuratWrappers); require(dplyr)
toLiger="species"
nmfK=30
args=commandArgs(trailingOnly = TRUE)
rdsObj=args[1]

runSeurat.umap.louvain=function(obj, ndim=10, louvainResolution=1.0 , reductionName="pca", tsne=FALSE){
  require(dplyr); require(Seurat)
  obj = obj %>% RunUMAP(reduction = reductionName, dims = 1:ndim) %>% 
    FindNeighbors(reduction = reductionName, dims = 1:ndim) %>% 
    FindClusters(resolution = louvainResolution)
  
  if (tsne){
    obj = obj %>% RunTSNE(reduction = reductionName, dims = 1:ndim)
  }
  return(obj)
}

reducedMNN_seurat=function(obj, nDim=10,  anchorFeature="batch",  clusterResolution=1.0, SEED=100, reduction="pca"){
  require(batchelor); require(Seurat); require(dplyr)
  set.seed(SEED)
  seuratPCs=obj@reductions[[reduction]]@cell.embeddings[,1:nDim]
  reducedMNN.obj = reducedMNN(seuratPCs, batch=obj@meta.data[[anchorFeature]])
  
  batchelor_reducedMNN_correctedPCs=reducedMNN.obj$corrected
  colnames(batchelor_reducedMNN_correctedPCs) = 1:ncol(batchelor_reducedMNN_correctedPCs) %>% paste0("reducedMNN_", .)
  obj[["reducedmnn"]] = CreateDimReducObject(embeddings = batchelor_reducedMNN_correctedPCs, key = "reducedMNN_", assay = DefaultAssay(obj))
  obj = obj %>% runSeurat.umap.louvain(., ndim=nDim, reductionName="reducedmnn", louvainResolution = clusterResolution)
  
  return(obj)
}

obj=readRDS(paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/",rdsObj))

obj = obj %>% NormalizeData(.) %>% FindVariableFeatures(.) %>%  
  ScaleData(., split.by = toLiger, do.center = FALSE) %>% 
  RunOptimizeALS(., k=nmfK, lambda=5, split.by = toLiger) %>%
  RunQuantileNorm(., split.by=toLiger)
obj = obj %>% RunUMAP(., dims = 1:ncol(.[["iNMF"]]), reduction = "iNMF")

obj_MNN = obj %>% reducedMNN_seurat(., nDim=nmfK, reduction="iNMF", anchorFeature = "batch")
saveRDS(obj, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/inmf_", rdsObj))
saveRDS(obj_MNN, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/inmf_mnn_", rdsObj))



