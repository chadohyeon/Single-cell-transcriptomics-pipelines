library(Seurat);library(infercnv); library(ggplot2); library(dplyr)

setwd("~/KHJ/RData/05.CNVanalysis/mouse_gbm2")

mousecnv <- subset(obj.final, idents = c("MG", "MAC", "EC", "PC", "pre-CC1", "pre-CC2", "CC1", "CC2", "CC3", "CC4", "NSC", "OPC", "COP", "OL", "TAC", "NB", "Neuron"))

table(mousecnv$orig.ident, mousecnv$celltype)

obj <- subset(mousecnv, downsample=200)


counts_matrix = obj@assays$RNA@counts[,colnames(obj)]
write.table(round(counts_matrix, digits=3), file='sc.10x.counts.matrix', quote=F, sep="\t")     

cellAnnotations=obj@meta.data %>% dplyr::select(celltype)
names(cellAnnotations) <- NULL
head(cellAnnotations)
write.table(cellAnnotations, file='m.cellAnnotations.txt', quote=F, sep="\t")      

mouse_gen_pos <- read.delim("/home/jungkim/KHJ/RData/05.CNVanalysis/mouse_gbm/mouse_gen_pos.txt", header = FALSE, row.names=1)
names(mouse_gen_pos) <- NULL
head(mouse_gen_pos)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=cellAnnotations,
                                    delim="\t",
                                    gene_order_file=mouse_gen_pos,
                                    ref_group_names=c("MG", "MAC", "EC", "PC"), chr_exclude=c('chrX', 'chrY'))

#new_gene_order = data.frame()
#for (chr_name in c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")) {
#  new_gene_order = rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])}
#names(new_gene_order) <- c("chr", "start", "stop")
#infercnv_obj@gene_order = new_gene_order
#infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)


######################################

setwd("~/KHJ/RData/05.CNVanalysis/mouse_gbm")
mousecnv2 <- subset(obj.final, idents = c("pre-CC1", "pre-CC2", "CC1", "CC2", "CC3", "CC4", "NSC", "OPC", "COP"))
Idents(mousecnv2) <- mousecnv2$orig.ident
table(mousecnv2$orig.ident, mousecnv2$celltype)

obj <- subset(mousecnv2, downsample=500)


counts_matrix = obj@assays$RNA@counts[,colnames(obj)]
write.table(round(counts_matrix, digits=3), file='sc.10x.counts.matrix', quote=F, sep="\t")     

cellAnnotations=obj@meta.data %>% dplyr::select(orig.ident)
names(cellAnnotations) <- NULL
head(cellAnnotations)
write.table(cellAnnotations, file='m.cellAnnotations.txt', quote=F, sep="\t")      

mouse_gen_pos <- read.delim("/home/jungkim/KHJ/RData/05.CNVanalysis/mouse_gbm/mouse_gen_pos.txt", header = FALSE, row.names=1)
names(mouse_gen_pos) <- NULL
head(mouse_gen_pos)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=cellAnnotations,
                                    delim="\t",
                                    gene_order_file=mouse_gen_pos,
                                    ref_group_names=c("Control"), chr_exclude=c('chrX', 'chrY'))



infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)














counts_matrix = GetAssayData(obj.final, slot="counts")
#counts_matrix = as.matrix(seurat_obj@assays$RNA@counts[,colnames(seurat_obj)])

# save the output table as an R object (faster and more size efficient)
saveRDS(round(counts_matrix, digits=3), "mouse.sc.10x.counts.matrix")
names(counts_matrix) <- NULL
# write the output table in txt format
write.table(round(counts_matrix, digits=3), file='sc.10x.counts.matrix', quote=F, sep="\t")         

cellAnnotations=obj.final@meta.data %>% dplyr::select(orig.ident)
names(cellAnnotations) <- NULL
head(cellAnnotations)
write.table(cellAnnotations, file='ms.cellAnnotations.txt', quote=F, sep="\t")      


### When using infercnv::run(), set 'cutoff=0.1' with 10xGenomics data, instead of the default (1) we tend to use with smartSeq2 and less sparse data.

cellAnnotations  <- read.delim("/home/jungkim/KHJ/RData/05.CNVanalysis/mouse_gbm/ms.cellAnnotations.txt", header = FALSE, row.names=1)
names(cellAnnotations) <- NULL
head(cellAnnotations)

mouse_gene_position <- read.delim("/home/jungkim/KHJ/RData/05.CNVanalysis/mouse_gbm
                                  /mouse_gen_pos.txt", header = FALSE, row.names=1)
#gene_ordering_file
#File > import dataset > txt
names(mouse_gene_position) <- NULL
head(mouse_gene_position)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=cellAnnotations,
                                    delim="\t",
                                    gene_order_file=mouse_gene_position ,
                                    ref_group_names=c("PC", "EC"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)
save(infercnv_obj, file="mouse_infercnv.RData")

