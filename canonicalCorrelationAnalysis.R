#!/usr/bin/Rscript
runSeuratPCA=function(obj, ndim=10, regress=c(), normalization.method="LogNormalize", sct_method="glmGamPoi"){
  require(dplyr); require(Seurat); require(glmGamPoi)
  if (normalization.method=="LogNormalize"){
    DefaultAssay(obj)="RNA"
    obj = obj %>% NormalizeData(.) %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% ScaleData(., vars.to.regress=regress)
  }else if(normalization.method=="SCT"){
    obj = obj %>% SCTransform(., method=sct_method, vars.to.regress=regress)
  }else if(normalization.method=="None"){
    obj = obj %>% ScaleData(., vars.to.regress=regress)
  }else{
    DefaultAssay(obj)="RNA"
    obj = obj %>% NormalizeData(.) %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% ScaleData(., vars.to.regress=regress)
  }
  obj = obj %>% RunPCA(., npcs=ndim)
}


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


wrapperRunHarmony=function(inputData,project_name="harmonyMerged", species="human", anchorFeature="orig.ident", ndim=10, louvainResolution=1.0, regress=c(), SCT_normalization=FALSE){
  require(harmony); require(dplyr); require(Seurat)
  if (mode(inputData)=="S4"){
    print("[Seurat S4] detected")
    data_vector=SplitObject(inputData, anchorFeature)
  }else if(mode(inputData)=="list" & unique(as.character(lapply(inputData, function(x) mode(x))))=="S4"){
    print("[list of Seurat S4] detected")
  }else{
    print("Please specify input data as [Seurat S4] or [list of Seurat S4s]")
    return(NULL)
  }
  cellIDs = data_vector %>% lapply(., function(x) x@project.name) %>% as.character()
  obj1=data_vector[[1]]
  obj2=data_vector[2:length(data_vector)]
  merged.obj = merge(obj1, y = obj2, add.cell.ids = cellIDs, project = project_name)
  merged.obj[["toAnchor"]] = merged.obj[[anchorFeature]]

  if(SCT_normalization){
    merged.obj = runSeuratPCA(merged.obj, ndim, regress, "SCT")
    merged.obj = merged.obj %>% RunHarmony("toAnchor", plot_convergence = TRUE, assay.use="SCT")
  }else{
    merged.obj = runSeuratPCA(merged.obj, ndim, regress, "LogNormalize")
    merged.obj = merged.obj %>% RunHarmony("toAnchor", plot_convergence = TRUE, assay.use="RNA")
  }
  
  merged.obj %>% DimPlot(object = ., reduction = "pca", group.by = "toAnchor", pt.size = .1) %>% print(.)
  merged.obj %>% VlnPlot(object = ., features = "PC_1", group.by = "toAnchor", pt.size = 0) %>% print(.)
  
  merged.obj %>% DimPlot(object = ., reduction = "harmony", group.by = "toAnchor", pt.size = .1) %>% print(.)
  merged.obj %>% VlnPlot(object = ., features = "harmony_1", group.by = "toAnchor", pt.size = 0) %>% print(.)
  
  # Downstream analysis
  merged.obj = runSeurat.umap.louvain(merged.obj, ndim, louvainResolution, reductionName="harmony")
  
  merged.obj %>% DimPlot(., reduction = "umap",label = TRUE, pt.size = .1) %>% print(.)
  merged.obj %>% DimPlot(., reduction = "umap",group.by = "toAnchor",label = TRUE, pt.size = .1) %>% print(.)
  return(merged.obj)
}

wrapperRunSeuratCCA=function(inputData,project_name="CCA_Merged", species="human", anchorFeature="orig.ident", fn="", ndim=30, louvainResolution=1.0, regress=c(), SCT_normalization=FALSE, k.weights=100, memMB=20000){
  require(dplyr); require(Seurat); require(glmGamPoi)
  options(future.globals.maxSize=memMB * 1024^2)
  if (mode(inputData)=="S4"){
    print("[Seurat S4] detected")
    data_vector=SplitObject(inputData, anchorFeature)
  }else if(mode(inputData)=="list" & unique(as.character(lapply(inputData, function(x) mode(x))))=="S4"){
    print("[list of Seurat S4] detected")
  }else{
    print("Please specify input data as [Seurat S4] or [list of Seurat S4s]")
    return(NULL)
  }
  cellIDs = data_vector %>% lapply(., function(x) x@project.name) %>% as.character()
  obj1=data_vector[[1]]
  obj2=data_vector[2:length(data_vector)]
  merged.obj = merge(obj1, y = obj2, add.cell.ids = cellIDs, project = project_name)

  individual_list = SplitObject(merged.obj, split.by = anchorFeature)
  minCellNum=individual_list %>% lapply(., function(x) length(Cells(x))) %>% unlist(.) %>% min(.)
  k.weights=min(k.weights, minCellNum)

  if(SCT_normalization){
    individual_list = lapply(X = individual_list, function(x) SCTransform(x, method="glmGamPoi", vars.to.regress = regress))
    SCT_features=SelectIntegrationFeatures(object.list = individual_list, nfeatures=3000) ### Seurat suggestion
    individual_list=PrepSCTIntegration(object.list = individual_list, anchor.features = SCT_features)
    CCA_anchors = FindIntegrationAnchors(object.list = individual_list, dims = 1:ndim, normalization.method = "SCT", anchor.features = SCT_features)
    saveRDS(CCA_anchors, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/anchor_cca_SCT_", fn))
    
    obj_integrated = IntegrateData(anchorset = CCA_anchors, dims = 1:ndim, normalization.method = "SCT", k.weight = k.weights)
    DefaultAssay(obj_integrated) = "integrated"
    obj_integrated = obj_integrated %>% RunPCA(., npcs=ndim)
  }else{
    individual_list = lapply(X = individual_list, FUN = NormalizeData)
    individual_list = lapply(X = individual_list, FUN = FindVariableFeatures)
    features = SelectIntegrationFeatures(object.list = individual_list)
    CCA_anchors = FindIntegrationAnchors(object.list = individual_list, dims = 1:ndim, anchor.features = features)
    saveRDS(CCA_anchors, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/anchor_cca_logNorm_", fn))
    obj_integrated = IntegrateData(anchorset = CCA_anchors, dims = 1:ndim, k.weight = k.weights)
    DefaultAssay(obj_integrated) = "integrated"
    obj_integrated = obj_integrated %>% runSeuratPCA(., ndim, regress, "None")
  }

  obj_integrated = obj_integrated %>% runSeurat.umap.louvain(., ndim=ndim, louvainResolution=louvainResolution)
  return(obj_integrated)
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


args=commandArgs(trailingOnly = TRUE)
rdsObj=args[1]
sctNorm=args[2]
k_weights=100

if(sctNorm=="y"){
  sctNorm=TRUE
  prep="SCT"
}else{
  sctNorm=FALSE
  prep="logNorm"}
merged_human_mouse=readRDS(paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/",rdsObj))

ccaRes = wrapperRunSeuratCCA(merged_human_mouse, species="human", anchorFeature = "species", SCT_normalization = sctNorm, k.weights=k_weights, fn=rdsObj)
ccaRes_MNN=ccaRes %>% reducedMNN_seurat(., nDim=30, reduction="pca", anchorFeature = "batch")

saveRDS(ccaRes, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/seurat_cca_",prep, "_", rdsObj))
saveRDS(ccaRes_MNN, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/seurat_cca_",prep, "_mnn_", rdsObj))
