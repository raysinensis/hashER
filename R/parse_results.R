parse_results <- function(res,
                          whitelist = NULL) {
  df <- data.table::fread(res)
  if (colnames(df)[1] == "V1") {
    # HashSolo
    df <- df %>% column_to_rownames("V1") %>% 
      rename(local = Classification) %>% 
      mutate(global = case_when(
        local == "Negative" ~ "Negative",
        local == "Doublet" ~ "Doublet",
        T ~ "Singlet"
      )) %>% 
      select(global, local)
  }
  if (colnames(df)[1] == "id") {
    # HTODemux
    df <- df %>% column_to_rownames("id") %>% 
      select(global, local)
  } else if (colnames(df)[1] == "barcodekey") {
    # demuxEM
    df <- df %>% rename(global = demux_type, local = assignment) %>% 
      mutate(global = str_to_title(global)) %>% 
      mutate(global = ifelse(global == "Unknown", "Negative", global)) %>% 
      mutate(local = case_when(
        global == "Doublet" ~ "Doublet",
        global == "Negative" ~ "Negative",
        TRUE ~ local)) %>% 
      column_to_rownames("barcodekey")
  }
  
  if (!is.null(whitelist)) {
    df <- df[whitelist, ]
  }
  
  df
}