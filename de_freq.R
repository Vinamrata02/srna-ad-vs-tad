nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

tbl <- pb_ds$table$TAD_vs_all

tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

top2 <- bind_rows(lapply(tbl_fil, top_n, 2, p_adj.loc))
format(top2[, -ncol(top2)], digits = 2)

frq <- calcExprFreqs(sce, assay = "counts", th = 0)

t(head(frq10))

gids <- levels(sce$group_id)
frq10 <- vapply(as.list(assays(frq)), 
                function(u) apply(u[, gids] > 0.1, 1, any), 
                logical(nrow(sce)))

tbl_fil2 <- lapply(kids, function(k)
  dplyr::filter(tbl_fil[[k]], 
                gene %in% names(which(frq10[, k]))))

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil2, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)