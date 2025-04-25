prepSCE <- function(x, 
                    kid = "ident", 
                    sid = "subject", 
                    gid = "group_id", 
                    drop = FALSE) {
  
  stopifnot(is(x, "SingleCellExperiment"))
  
  args <- as.list(environment())
  ids <- args[grep("[a-z]id", names(args))]
  ids <- unlist(ids)
  
  stopifnot(is.character(ids))
  stopifnot(all(ids %in% colnames(colData(x))))
  
  cd0 <- colData(x)
  cd <- data.frame(cd0[ids], check.names = FALSE)
  cd <- mutate_all(cd, as.factor)
  colnames(cd) <- unlist(formals()[names(ids)])
  
  if (!drop)
    cd <- data.frame(cd,
                     cd0[setdiff(colnames(cd0), ids)], 
                     check.names = FALSE)
  
  # replace colData in SCE
  colData(x) <- DataFrame(cd)
  
  # construct metadata
  ei <- .make_ei(x)
  metadata(x)$experiment_info <- ei
  
  return(x)
}