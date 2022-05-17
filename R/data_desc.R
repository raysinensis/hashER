hash_data <- setNames(
  data.frame(matrix(ncol = 4, nrow = 0)), 
  c("sample_id","id_col", "ids", "code")
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "A",
    id_col = "sex_dmg",
    ids = "F_M",
    code = "str_sub(id_hash, 1, 1)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "P",
    id_col = "freemuxlet_id",
    ids = "3911_9319",
    code = "str_sub(id_hash, 3, 6)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "hash18",
    id_col = "type",
    ids = "HEK_K562_KG1_THP1",
    code = "str_remove(id_hash, '[\\.-].+')"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "2hash",
    id_col = "type",
    ids = "line_pbmc",
    code = ""
  )
)