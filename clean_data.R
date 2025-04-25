library(tidyverse)
library(muscat)
library(Seurat)
library(limma)
library(SingleCellExperiment)

# read clusterwise (inhibitory) annotated data (Seurat object)
object = readRDS("/N/slate/vinshar/Excitatory_neurons_set1.rds")

# read label data shared by Prof (R data frame)
labels <- read_csv("/N/u/vinshar/Quartz/Downloads/ROSMAP_Clinical_Subtypes(in)(1).csv")

# keep only the subject and Subtype column
labels <- labels[, c("subject", "Subtype")]
# keep only TAD and AsymAD
labels <- labels[labels$Subtype %in% c("TAD", "AsymAD"), ]

metadata <- object@meta.data
metadata <- metadata[, c("subject", "cell_type_high_resolution")]
rownames <- rownames(object@meta.data)
metadata <- left_join(metadata, labels)
rownames(metadata) <- rownames
object@meta.data <- metadata 
rows_with_na <- rownames(object@meta.data)[apply(object@meta.data, 1, function(x) any(is.na(x)))]
object@assays$RNA@counts <- object@assays$RNA@counts[, !colnames(object@assays$RNA@counts) %in% rows_with_na]
object@assays$RNA@data <- object@assays$RNA@data[, !colnames(object@assays$RNA@data) %in% rows_with_na]
object@reductions <- list()
object@meta.data <- na.omit(object@meta.data)
object <- SetIdent(object, value=object@meta.data$cell_type_high_resolution)

# 1% subset
# sample_size <- ceiling(0.2 * ncol(object))
# selected_cells <- sample(colnames(object), size=sample_size, replace=FALSE)
# object_sample <- subset(object, cells=selected_cells)

# convert Seurat object to SCE because muscat operates on SCE objects
sce <- as.SingleCellExperiment(object)

# subset
#sce <- sce[, sce$group_id %in% c("TAD", "AsymAD")]

# rename columns according to what muscat expects
sce <- prepSCE(sce, kid="cell_type_high_resolution", sid="subject", gid="Subtype")

# pseudo-bulk analysis
pb <- aggregateData(sce, assay="counts", fun="sum", by=c("cluster_id", "sample_id"))

# create contrast groups
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
# contrast <- makeContrasts(TAD_vs_CN = TAD - CN, TAD_vs_LowNFT = TAD - LowNFT, TAD_vs_AsymAD = TAD - AsymAD, CN_vs_LowNFT = CN - LowNFT, CN_vs_AsymAD = CN - AsymAD, LowNFT_vs_AsymAD = LowNFT - AsymAD, levels=mm)

# one-vs-all
# contrast <- makeContrasts(TAD_vs_all = 3*TAD - (CN + LowNFT + AsymAD),
                          #CN_vs_all = 3*CN - (TAD + LowNFT + AsymAD),
                          #LowNFT_vs_all = 3*LowNFT - (TAD + CN + AsymAD),
                          #AsymAD_vs_all = 3*AsymAD - (TAD + CN + LowNFT),
                          #levels = mm)

# TAD vs AsymAD
contrast <- makeContrasts(TAD_vs_AsymAD = TAD - AsymAD, levels=mm)

# test
# contrast <- makeContrasts(TAD_vs_CN = TAD - CN, TAD_vs_LowNFT = TAD - LowNFT, TAD_vs_AsymAD = TAD - AsymAD, CN_vs_LowNFT = CN - LowNFT, CN_vs_AsymAD = CN - AsymAD, LowNFT_vs_AsymAD = LowNFT - AsymAD, levels=mm)


# DE
pb_ds <- pbDS(pb, method="edgeR", design = mm, contrast = contrast)
#pb_ds <- pbDS(pb, method="edgeR")

# Plot
#pbHeatmap(sce, pb_ds, top_n=5)
pbHeatmap(sce, pb_ds, fdr=0.05, lfc=1, top_n=5, c="TAD_vs_AsymAD", sort_by = "logFC", decreasing = TRUE)

# analysis
# count
vapply(pb_ds$table$TAD_vs_AsymAD, function(u) sum (abs(u$logFC) > 1 & u$p_adj.loc < 0.05),  numeric(1))
# genes
lapply(pb_ds$table$TAD_vs_AsymAD, function(u) {u %>% filter((abs(logFC) > 1) & p_adj.loc < 0.05) %>% arrange(p_adj.loc)})
                                                                                                      
#to convert the DEG in csv file
csv <- lapply(pb_ds$table$TAD_vs_AsymAD, function(u) {u %>% filter((abs(logFC) > 1) & p_adj.loc < 0.05) %>% arrange(p_adj.loc)})
                                                                                                      
df_list <- lapply(csv, function(x) {
  if (nrow(x) > 0) return(as.data.frame(x))  # Keep only non-empty elements
})
                                                                                                      
# Bind all non-empty data frames into one
df_combined <- do.call(rbind, df_list)
                                                                                                      
# write
write.csv(df_combined, "/N/u/vinshar/Quartz/TAD_ASYMAD/csv_TAD_ASYMAD_exc1.csv")