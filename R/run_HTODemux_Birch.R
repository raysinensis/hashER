run_HTODemux_Birch <- function(raw_counts = NULL,
                         filtered_counts = NULL,
                         whitelist = NULL,
                         whitelist_rm = "-1$",
                         threshold = 100,
                         norm_total = 10000,
                         kfunc = "birch",
                         out_name = "HTODemux.csv.gz") {
  # get counts
  if (!is.null(whitelist)) {
    whitelist <- readr::read_lines(whitelist)
    if (!is.null(whitelist_rm)) {
      whitelist <- whitelist %>% str_remove(whitelist_rm)
    }
    message("using ", length(whitelist), " cell ids")
  }
  if (!is.null(raw_counts)) {
    
    counts <- data.table::fread(raw_counts) %>% column_to_rownames("V1")
  } else if (!is.null(filtered_counts)) {
    counts <- data.table::fread(filtered_counts) %>% column_to_rownames("V1")
  } else {
    stop("no counts specified")
  }
  
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
  counts_norm <- sweep(counts,1,rowSums(counts),FUN="/")
  counts_norm <- round(counts_norm*1000)
  so <- Seurat::CreateSeuratObject(counts = t(as.matrix(counts_norm)), assay="HTO")
  so <- Seurat::NormalizeData(so, normalization.method = "CLR")
  so <- HTODemux_mod(so, kfunc = kfunc)
  
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

HTODemux_mod <- function(
  object,
  assay = "HTO",
  positive.quantile = 0.99,
  init = NULL,
  nstarts = 100,
  kfunc = "clara",
  nsamples = 100,
  seed = 42,
  verbose = TRUE,
  niter = 10
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  #initial clustering
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(
    object = object,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(
    EXPR = kfunc,
    'kmeans' = {
      init.clusters <- kmeans(
        x = t(x = GetAssayData(object = object, assay = assay)),
        centers = ncenters,
        nstart = nstarts
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    },
    'clara' = {
      #use fast k-medoid clustering
      init.clusters <- cluster::clara(
        x = t(x = GetAssayData(object = object, assay = assay)),
        k = ncenters,
        samples = nsamples
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
    },
    'hdbscan' = {
      init.clusters <- dbscan::hdbscan(
        x = t(x = GetAssayData(object = object, assay = assay)),
        minPts = nsamples,
        verbose = T,
        gen_simplified_tree = T
      )
      #identify positive and negative signals for all HTO
      res <<- init.clusters
      print(res)
      Idents(object = object, cells = names(x = init.clusters$cluster), drop = TRUE) <- init.clusters$cluster
    },
    'birch' = {
      i = 0.5
      message("tuning...")
      
      while (i > 0) {
        dat <- t(x = GetAssayData(object = object, assay = assay))
        birch <- stream::DSC_BIRCH(threshold = i, branching = 5, maxLeaf = 5, outlierThreshold = 2*i)
        ds <- stream::DSD_Memory(dat, loop = TRUE)
        update(birch, ds, n = 3 * nrow(dat), verbose = T, block = 3 * nrow(dat))
        init.clusters <- stream::get_assignment(birch, stream::get_points(ds, n = nrow(dat)))
        if (length(unique(init.clusters)) <= ncol(dat)) {
          i <- i - 0.005
        } else if (length(unique(init.clusters)) > 10 * ncol(dat)) {
          i <- i + 0.21
        } else if ((length(unique(init.clusters)) > 1.1 * ncol(dat)) & (length(unique(init.clusters)) >  1 + ncol(dat))) {
          i <- i + 0.032
        } else {
          break
        }
      }
      print(table(init.clusters))
      Idents(object = object, cells = names(x = init.clusters), drop = TRUE) <- init.clusters
    },
    'gmm' = {
      dat <- t(x = GetAssayData(object = object, assay = assay))
      gmm <- ClusterR::GMM(dat, 
                           gaussian_comps = ncol(dat) + 1, 
                           dist_mode = "maha_dist", 
                           seed_mode = "static_spread",
                           km_iter = niter,
                           em_iter = niter,
                           seed = seed
      )
      init.clusters <- predict(gmm, newdata = dat)
      print(table(init.clusters))
      Idents(object = object, cells = names(x = init.clusters), drop = TRUE) <- init.clusters
    },
    stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
  )
  #average hto signals per cluster
  average.expression <- AverageExpression(
    object = object,
    assays = assay,
    verbose = FALSE
  )[[assay]]
  
  # average.expression <- clustifyr::average_clusters(
  #   object[[assay]]@data,
  #   metadata = Idents(object),
  #   if_log = TRUE,
  #   output_log = FALSE,
  #   method = "trimean"
  # )
  # 
  #checking for any cluster with all zero counts for any barcode
  if (sum(average.expression == 0) > 0) {
    warning("Cells with zero counts exist as a cluster.")
  }
  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    values.use <- values[WhichCells(
      object = object,
      idents = levels(x = Idents(object = object))[sort(average.expression[iter, ], index.return = TRUE)$ix[1:max(round(length(average.expression[iter, ])*0.2), 1)]]
      # idents = "0"
      # idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, ])]]
    )]
    fit <- suppressWarnings(expr = fitdistrplus::fitdist(data = values.use, distr = "nbinom", lower = c(0,0)))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, " reads"))
    }
  }
  # now assign cells to HTO based on discretized values
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = Seurat:::MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.max[x])[1])
    }
  )])
  hash.secondID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.second[x])[1])
    }
  )])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(
    X = 1:length(x = hash.maxID),
    FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
    }
  )
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(
    hash.maxID,
    hash.secondID,
    hash.margin,
    classification,
    classification.global
  )
  colnames(x = classification.metadata) <- paste(
    assay,
    c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
    sep = '_'
  )
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, '_classification')
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- 'Doublet'
  object$hash.ID <- Idents(object = object)
  return(object)
}
