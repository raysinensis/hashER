run_HTODemux <- function(raw_counts = NULL,
                         filtered_counts = NULL,
                         whitelist = NULL,
                         whitelist_rm = "-1$",
                         threshold = 100,
                         out_name = "HTODemux.csv.gz",
                         use_mod = F) {
  # get counts
  if (!is.null(whitelist)) {
    whitelist <- readr::read_lines(whitelist)
    if (!is.null(whitelist_rm)) {
      whitelist <- whitelist %>% str_remove(whitelist_rm)
    }
    message("using ", length(whitelist), " cell ids")
  }
  if (!is.null(raw_counts)) {
    counts <- data.table::fread(raw_counts)
    counts <- counts %>% column_to_rownames(colnames(counts)[1])
  } else if (!is.null(filtered_counts)) {
    counts <- data.table::fread(filtered_counts)
    counts <- counts %>% column_to_rownames(colnames(counts)[1])
  } else {
    stop("no counts specified")
  }

  message("HTO dims: ", paste0(dim(counts), collapse = ","))
  # fill in missing cells
  if (!is.null(whitelist)) {
    emptys <- setdiff(whitelist, rownames(counts))
    if (length(emptys) > 0) {
      empty_mat <- matrix(data = 0, ncol = ncol(counts), nrow = length(emptys))
      rownames(empty_mat) <- emptys
      colnames(empty_mat) <- colnames(counts)
      counts <- rbind(counts, empty_mat)
    }
    counts <- counts[whitelist, ]
  }

    # remove low count drops first
  negative_by_default <- rownames(counts)[rowSums(counts) < threshold]
  counts <- counts[rowSums(counts) >= threshold,]
  so <- Seurat::CreateSeuratObject(counts = t(as.matrix(counts)), assay="HTO")
  so <- Seurat::NormalizeData(so, normalization.method = "CLR")
  so <- Seurat::HTODemux(so)
  
  # results out
  HTO_results <- so@meta.data %>% 
    rownames_to_column("id") %>% 
    select(id, global = HTO_classification.global, local = hash.ID, 
           HTO_classification, HTO_maxID, HTO_secondID, HTO_margin)
  if (length(negative_by_default) > 0) {
    HTO_neg <- data.frame(id = negative_by_default, 
                          global = "Negative", 
                          local = "Negative", 
                          HTO_classification = "Negative",
                          HTO_maxID = "low", HTO_secondID = "low", HTO_margin = 0)
    
    HTO_results <- bind_rows(HTO_results, HTO_neg) %>% 
      column_to_rownames("id")
  } else {
    HTO_results <- HTO_results %>% 
      column_to_rownames("id")
  }
  
  # write results
  if (!is.null(whitelist)) {
    HTO_results <- HTO_results[whitelist, ]
  }
  HTO_results %>% rownames_to_column("id") %>% 
    readr::write_csv(out_name)
  HTO_results
}
