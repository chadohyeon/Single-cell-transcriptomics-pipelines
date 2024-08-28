#!/usr/bin/Rscript

require(magrittr); require(tradeSeq); require(slingshot); require(dplyr); require(monocle3)

args=commandArgs(trailingOnly = TRUE)
fn=args[1]
#path="/home/dcha/02.glioblastoma_scRNAseq/rdata/monocle3_data/"
path="/home/jungkim/KHJ/RData/02.subset_renewed/obj1_monocle3.rds"
cds=readRDS(paste0(path,fn))

cds <- obj2_monocle3
mst = principal_graph(cds)$UMAP

y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

# Get the root vertices
# It is the same node as above
root = cds@principal_graph_aux$UMAP$root_pr_nodes

# Get the other endpoints
endpoints = names(which(igraph::degree(mst) == 1))
endpoints = endpoints[!endpoints %in% root]

# For each endpoint
cellWeights = lapply(endpoints, function(endpoint) {
  # We find the path between the endpoint and the root
  path = igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path = as.character(path)
  # We find the cells that map along that path
  df = y_to_cells[y_to_cells$Y %in% path, ]
  df = data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) = endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()
rownames(cellWeights) = colnames(cds)
pseudotime = matrix(pseudotime(cds), ncol = ncol(cellWeights),
                    nrow = ncol(cds), byrow = FALSE)
#icMat = evaluateK(counts = cds@assays@data$counts, pseudotime = pseudotime, cellWeights = cellWeights,
#                   k=3:7, nGenes = 100, verbose = TRUE, plot = TRUE)
fitGAM_cds = fitGAM(counts = cds@assays@data$counts,
                    pseudotime = pseudotime,
                    cellWeights = cellWeights)

fitGAM_cds_fn = fn %>% gsub(".rds", ".fitGAM.rds", .)
saveRDS(fitGAM_cds, paste0(path, fitGAM_cds_fn))

saveRDS(fitGAM_cds, file="obj2_fitGAM.rds")




########################## 
#https://bioconductor.org/packages/devel/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html#within-lineage-comparisons
#Within-lineage comparisons

fitGAM_cds <- obj2_fitGAM
assoRes <- associationTest(fitGAM_cds)
head(assoRes)

startRes <- startVsEndTest(fitGAM_cds)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(fitGAM_cds)[oStart[3]]


countsED<-cds@assays@data@listData$counts
plotSmoothers(fitGAM_cds, countsED, gene = sigGeneStart)

plotGeneCount(crv, counts, gene = sigGene)


