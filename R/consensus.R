find_consensus <- function(df,
                           cols = contains("_local"),
                           new_col = "consensus") {
  df2 <- df %>% select(all_of(cols))
  df2[[paste0(new_col, "_local")]] <- apply(df2, MARGIN = 1, FUN = function(x) find_singlet(x))
  df2[[paste0(new_col, "_global")]] <- ifelse(df2[[paste0(new_col, "_local")]] %in% c("Negative", "Doublet", "Discord"), df2[[paste0(new_col, "_local")]], "Singlet")
  cbind(df, df2[, c(paste0(new_col, "_local"), paste0(new_col, "_global"))])
}

find_singlet <- function(vec, discord = "Discord") {
  vec <- unique(vec)
  if (length(vec) == 1) {
    return(vec)
  }
  vec <- vec[!(vec == "Negative")]
  if (length(vec) == 1) {
    return(vec[1])
  }
  vec <- vec[!(vec == "Doublet")]
  if (length(vec) == 1) {
    return(vec[1])
  } else {
    return(discord)
  }
}