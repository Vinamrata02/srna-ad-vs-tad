library(UpSetR)
library(scater)

# upsetset plot
tbl <- pb_ds$table[[1]]
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
  })

de_gs_by_k <- map(tbl_fil, "gene")
upset(fromList(de_gs_by_k))

# upset plot from merged
df <- read.csv('/N/slate/vinshar/AsymAD_vs_ALL_merged_des.csv')
tbl <- split(df, df$cell_type)
de_gs_by_k <- map(tbl, "gene")
upset(fromList(de_gs_by_k))

