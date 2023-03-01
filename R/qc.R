# return pvalue of bimodal test
test_bimodal <- function(mat,
                         threshold = 1,
                         return = "summary") {
  mat <- mat[rowSums(mat) >= threshold, ]
  if (nrow(mat) > 5000) {
    mat <- mat[sample(1:nrow(mat), 5000), ]
  }
  res_lmoc <- apply(
    mat,
    MARGIN = 2, 
    FUN = function(x) {
      suppressWarnings(multimode::modetest(x, mod0 = 1))
    }
  )
  
  if (return == "all") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$p.value
    }) %>% unlist()
    return(list(res_lmoc, mean(res_lmoc2)))
  } else if (return == "summary") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$p.value
    }) %>% unlist()
    return(mean(res_lmoc2))
  }
}

# return max_background_tag/min_background_tag ratio
test_background <- function(mat,
                            threshold = 1,
                            q = 0.5,
                            return = "summary") {
  mat <- mat[rowSums(mat) >= threshold, ]
  mat2 <- apply(
    mat, 
    MARGIN = 2,
    FUN = function(x) {
      x[x <= quantile(x, q)]
    }
  )
  lmin <- map(mat2, length) %>% unlist() %>% min()
  mat2 <- map(mat2, function(x) {x[1:lmin]})
  mat2 <- do.call(cbind, mat2)
  res_lmoc <- apply(
    mat2,
    MARGIN = 2, 
    FUN = function(x) {
      if (max(x) > 0) {
        suppressWarnings(multimode::locmodes(x))
      } else {
        list(locations = 0.001)
      }
    }
  )
  
  if (return == "all") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$locations
    }) %>% unlist()
    res_lmoc2 <- max(res_lmoc2) / min(res_lmoc2)
    return(list(res_lmoc, res_lmoc2))
  } else if (return == "summary") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$locations
    }) %>% unlist()
    return(max(res_lmoc2) / min(res_lmoc2))
  }
}

test_typebias <- function(mat,
                          types,
                          return = "summary") {
  # mat <- mat[rowSums(mat) >= threshold, ]
  cells <- split(rownames(mat), types)
  res <- map(cells, function(x) {
    temp <- mat[x,]
    apply(
      temp,
      MARGIN = 2, 
      FUN = function(x) {
        suppressWarnings(multimode::locmodes(x, mod0 = 2))
      }
    )
  })

  if (return == "all") {
    res2 <- map(res, function(y) {
      map(y, function(x) {
        x$locations[2]
      })
    } %>% unlist()) %>% as.data.frame() %>% t()
    res2 <- sum((matrixStats::colMaxs(res2) / matrixStats::colMins(res2)) >= 5) / ncol(res2)
    return(list(res, res2))
  } else if (return == "summary") {
    res2 <- map(res, function(y) {
      map(y, function(x) {
        x$locations
      })
    } %>% unlist()) %>% as.data.frame() %>% t()
    return(sum((matrixStats::colMaxs(res2) / matrixStats::colMins(res2)) >= 5) / ncol(res2))
  }
}

test_lowreads <- function(mat,
                     threshold = 1,
                     q = 0.99,
                     return = "summary") {
  mat <- mat[rowSums(mat) >= threshold, ]
  mat2 <- apply(
    mat, 
    MARGIN = 2,
    FUN = function(x) {
      x[x <= quantile(x, q)]
    },
    simplify = F
  )
  lmin <- map(mat2, length) %>% unlist() %>% min()
  mat2 <- map(mat2, function(x) {x[1:lmin]})
  mat2 <- do.call(cbind, mat2)
  res_lmoc <- apply(
    mat2,
    MARGIN = 2, 
    FUN = function(x) {
      suppressWarnings(multimode::locmodes(x, mod0 = 2))
    }
  )
  
  if (return == "all") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$locations[2]
    }) %>% unlist()
    return(list(res_lmoc, mean(res_lmoc2 <= 50)))
  } else if (return == "summary") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$locations
    }) %>% unlist()
    return(mean(res_lmoc2 <= 50))
  }
}

get_saturation <- function(mat,
                           logfile) {
  nmol <- sum(mat)
  if (str_to_lower(strsplit(basename(logfile), split="\\.")[[1]][-1]) == "json") {
    # parse kite
    nmap <- rjson::fromJSON(file = logfile)$n_pseudoaligned
  } else if (str_to_lower(strsplit(basename(logfile), split="\\.")[[1]][-1]) == "csv") {
    # parse cellranger
    # old version
    
    nmap <- tryCatch(
      read_csv(logfile) %>%
        filter(`Metric Name` == "Number of reads",
               `Library Type` == "Multiplexing Capture", 
               `Grouped By` == "Physical library ID") %>%
        pull(`Metric Value`) %>% 
        str_remove_all(","),
      error = function(e) {
        read_csv(logfile) %>%
          pull(`Antibody: Number of Reads`)
      }
    )
  } else if (str_to_lower(strsplit(basename(logfile), split="\\.")[[1]][-1]) == c("fastq", "gz")) {
    # parse fastq.gz
    nlines <- 0
    for (l1 in logfile) {
      f <- gzfile(l1, open="rb")
      while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
        nlines <- nlines + sum(chunk == as.raw(10L))
      }
    }
    nmap <- nlines/4
  }
  return(1 - nmol/as.numeric(nmap))
}

read_mat <- function(raw_counts = NULL,
                        filtered_counts = NULL,
                        whitelist = NULL,
                        whitelist_rm = "-1$") {
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
  
  counts
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x+1), na.rm=na.rm) / length(x)) - 1
}