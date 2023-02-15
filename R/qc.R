test_bimodal <- function(mat,
                         threshold = 10,
                         return = "summary") {
  mat <- mat[rowSums(mat) >= threshold, ]
  if (nrow(mat) > 20000) {
    mat <- mat[sample(1:nrow(mat), 20000), ]
  }
  res_lmoc <- apply(
    mat,
    MARGIN = 2, 
    FUN = function(x) {
      multimode::modetest(x)
    }
  )
  
  if (return == "all") {
    return(res_lmoc)
  } else if (return == "summary") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$p.value
    }) %>% unlist()
    return(mean(res_lmoc2))
  }
}

test_background <- function(mat,
                            threshold = 10,
                            q = 0.1,
                            return = "summary") {
  mat <- mat[rowSums(mat) >= threshold, ]
  res_lmoc <- apply(
    mat,
    MARGIN = 2, 
    FUN = function(x) {
      multimode::locmodes(x)
    }
  )
  
  if (return == "all") {
    return(res_lmoc)
  } else if (return == "summary") {
    res_lmoc2 <- map(res_lmoc, function(x) {
      x$locations
    }) %>% unlist()
    return(max(res_lmoc2) / min(res_lmoc2) >= 10)
  }
}

test_typebias <- function(mat,
                          types,
                          threshold = 10,
                          return = "summary") {
  mat <- mat[rowSums(mat) >= threshold, ]
  cells <- split(Cells(so), so$type)
  res <- map(cells, function(x) {
    temp <- mat[x,]
    apply(
      temp,
      MARGIN = 2, 
      FUN = function(x) {
        multimode::locmodes(x, mod0 = 2)
      }
    )
  })

  if (return == "all") {
    return(res)
  } else if (return == "summary") {
    res2 <- map(res, function(y) {
      map(y, function(x) {
        x$locations[2]
      })
    } %>% unlist()) %>% as.data.frame() %>% t()
    return(sum(matrixStats::colMaxs(res2) / matrixStats::colMins(res2) >=5) / ncol(res2))
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