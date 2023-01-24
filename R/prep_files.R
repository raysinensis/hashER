prep_files <- function(hto_mtx,
                       hto_feature,
                       hto_bc,
                       rna_mtx,
                       rna_feature,
                       rna_bc,
                       dir) {
  m <- Matrix::readMM(hto_mtx)
  m <- t(as.matrix(m))
  colnames(m) <- read_tsv(hto_feature, col_names = F) %>% pull("X1")
  rownames(m) <- read_tsv(hto_bc, col_names = F) %>% pull("X1")
  m2 <- Matrix::readMM(rna_mtx)
  rownames(m2) <- read_tsv(rna_feature, col_names = F) %>% pull("X1")
  colnames(m2) <- read_tsv(rna_bc, col_names = F) %>% pull("X1") %>% str_remove("-1")
  so <- CreateSeuratObject(m2)
  so <- so %>% NormalizeData() %>% 
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures        = 2000
    ) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(verbose = FALSE, resolution = 0.1) %>%
    RunUMAP(reduction = "pca", dims = 1:30)
  dir.create(dir,
             showWarnings = FALSE, 
             recursive = TRUE)
  message("Saving Seurat object...")
  saveRDS(so, file.path(dir, "so.rds"))
  message("Saving whitelist...")
  bcs <- intersect(colnames(m2), rownames(m))
  write_lines(bcs, file.path(dir, "barcodes.tsv.gz"))
  message("Saving HTO counts table...")
  data.table::fwrite(as.data.frame(m),
                     file.path(dir, "HTO_counts.csv.gz"),
                     row.names = TRUE, quote = FALSE)
  message("Done.")
}