#!/usr/bin/Rscript
require(dplyr); require(Seurat)

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

cca_anchor_transfer=function(query_obj, ref_obj, nDim=30, mnn_feature="batch"){
  require(Seurat); require(dplyr)
  ref_obj=ref_obj %>% runSeuratPCA(., ndim=nDim, normalization.method = "SCT")
  query_obj=query_obj %>% runSeuratPCA(., ndim=nDim, normalization.method = "SCT")
  query_obj[["pca"]]=NULL
  
  musHom.anchors = FindTransferAnchors(reference = ref_obj, query = query_obj, features = VariableFeatures(object = ref_obj),
                                       reference.assay = "SCT", query.assay = "SCT", reduction = "cca", normalization.method = "SCT")
  predictions = TransferData(anchorset = musHom.anchors, refdata = ref_obj$celltype,
                             dims = 1:nDim, weight.reduction="cca")
  query_obj$celltype_prev=query_obj$celltype
  query_obj$celltype=predictions[["predicted.id"]]
  
  
  #ref_obj_mnn = ref_obj %>% reducedMNN_seurat(., nDim=nDim, anchorFeature=mnn_feature)
  #ref_obj_mnn = ref_obj_mnn %>% RunUMAP(., dims=1:nDim, reduction="reducedmnn", return.model = TRUE)
  #query_obj_mnn = MapQuery(anchorset=musHom.anchors, reference=ref_obj_mnn, query=query_obj,
  #                     refdata = list(celltype = "celltype"), reference.reduction = "reducedmnn", reduction.model = "umap") 
  
  ref_obj = ref_obj %>% RunUMAP(., dims=1:nDim, reduction="pca", return.model = TRUE)
  query_obj = MapQuery(anchorset=musHom.anchors, reference=ref_obj, query=query_obj,
                       refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap") 
  
  resList=list(query_obj, ref_obj)
  names(resList)=c("query", "ref")  
  #resList=list(query_obj, ref_obj, query_obj_mnn, ref_obj_mnn)
  #names(resList)=c("query", "ref", "query_mnn", "ref_mnn")
  return(resList)
}

reducedMNN_seurat=function(obj, nDim=10,  anchorFeature="batch",  clusterResolution=1.0, SEED=100, reduction="pca"){
  require(batchelor); require(Seurat); require(dplyr)
  set.seed(SEED)
  seuratPCs=obj@reductions[[reduction]]@cell.embeddings[,1:nDim]
  reducedMNN.obj = reducedMNN(seuratPCs, batch=obj@meta.data[[anchorFeature]])
  
  batchelor_reducedMNN_correctedPCs=reducedMNN.obj$corrected
  colnames(batchelor_reducedMNN_correctedPCs) = 1:ncol(batchelor_reducedMNN_correctedPCs) %>% paste0("reducedMNN_", .)
  obj[["reducedmnn"]] = CreateDimReducObject(embeddings = batchelor_reducedMNN_correctedPCs, key = "reducedMNN_", assay = DefaultAssay(obj))

  return(obj)
}

args=commandArgs(trailingOnly = TRUE)
rdsObj=args[1]
refName="mouse"
queryName="human"

obj=readRDS(paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/musHom/",rdsObj))

obj.list=SplitObject(obj, "species")
ref_obj = obj.list[[refName]]
query_obj = obj.list[[queryName]]
remove(obj)
remove(obj.list)

resList = cca_anchor_transfer(query_obj = query_obj, ref_obj = ref_obj)

saveRDS(resList[["query"]], paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cca_transfer/cca_sct_Trans_",queryName,"_", rdsObj))
saveRDS(resList[["ref"]], paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cca_transfer/cca_sct_Trans_",refName,"_", rdsObj))
#saveRDS(resList[["query_mnn"]], paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cca_transfer/cca_sct_Trans_mnn_",queryName,"_", rdsObj))
#saveRDS(resList[["ref_mnn"]], paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cca_transfer/cca_sct_Trans_mnn_",refName,"_", rdsObj))
