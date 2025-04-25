library(scater)

tbl <- pb_ds$table$TAD_vs_all

tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

plotExpression(sce[, sce$cluster_id == "Inh RYR3 TSHZ2"],
               features = tbl_fil$`Inh RYR3 TSHZ2`$gene[seq_len(6)],
               x = "sample_id", colour_by = "group_id", ncol = 3) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

