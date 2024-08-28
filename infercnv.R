library(Seurat);library(infercnv); library(ggplot2); library(dplyr)


############################## human
#DefaultAssay(Hom.query) <- "RNA"
Idents(Hom.query) <- Hom.query$predicted.id
table(Hom.query$batch, Hom.query$predicted.id)
#humancnv <- subset(Hom.query, idents=c("pre-CC1", "pre-CC2", "CC1", "CC2", "CC3", "CC4", "NSC", "TAC", "NB", "Neuron", "OPC", "COP", "OL", "EPC", "EC", "PC", "VSMC", "MG", "MAC", "Meninges")) 
humancnv <- subset(Hom.query, idents=c("TAM1", "TAM2"), invert = TRUE) 

humancnv@active.ident <- factor(x = humancnv@active.ident, 
                                  levels = c("pre-CC1", "pre-CC2", "CC1", "CC2", "CC3", "CC4", "NSC", "TAC", "NB", "Neuron", "OPC", "COP", "OL", "EPC", "EC", "PC", "VSMC", "MG", "MAC", "Meninges"))

humancnv$celltype2 <- Idents(humancnv)

table(humancnv$batch, humancnv$celltype2)

Idents(humancnv) <- humancnv$region

humansvz <- subset(humancnv, idents="SVZ")
humantu <- subset(humancnv, idents="Tumor")

Idents(humansvz) <- humansvz$celltype2

table(humansvz$celltype2)
#humansvz <- subset(humansvz, downsample=200)


############
Idents(humantu) <- humantu$celltype2

table(humantu$celltype2)
humantu <- subset(humantu, downsample=200)



obj <- humansvz

counts_matrix = obj@assays$RNA@counts[,colnames(obj)]


# save the output table as an R object (faster and more size efficient)
#saveRDS(round(counts_matrix, digits=3), "sc.10x.counts.matrix")

# write the output table in txt format
write.table(round(counts_matrix, digits=3), file='sc.10x.counts.matrix', quote=F, sep="\t")     



cellAnnotations=obj@meta.data %>% dplyr::select(celltype2)
names(cellAnnotations) <- NULL
head(cellAnnotations)
write.table(cellAnnotations, file='h.cellAnnotations.txt', quote=F, sep="\t")      


#gene_ordering_file (File > import dataset > txt)

human_gen_pos <- read.delim("/home/jungkim/KHJ/RData/08.CNVanalysis/human_gen_pos.txt", header = FALSE, row.names=1)
names(human_gen_pos) <- NULL
head(human_gen_pos)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=cellAnnotations,
                                    delim="\t",
                                    gene_order_file=human_gen_pos,
                                    ref_group_names=c("MG", "MAC", "OL"), chr_exclude=c('X', 'Y', 'MT'))

#https://github.com/broadinstitute/infercnv/issues/290
new_gene_order = data.frame()
for (chr_name in c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")) {
  new_gene_order = rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
infercnv_obj@gene_order = new_gene_order
infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]


# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T, resume_mode=FALSE)


