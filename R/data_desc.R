hash_data <- setNames(
  data.frame(matrix(ncol = 4, nrow = 0)), 
  c("sample","id_col", "ids", "code")
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "A17k",
    id_col = "sex_dmg",
    ids = "F_M",
    code = "str_sub(id_hash, 1, 1)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "A45k",
    id_col = "sex_dmg",
    ids = "F_M",
    code = "str_sub(id_hash, 1, 1)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "P1mc",
    id_col = "freemuxlet_id",
    ids = "3911_9319",
    code = "str_sub(id_hash, 3, 6)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "P2a",
    ids = "3911_9319",
    id_col = "freemuxlet_id",
    code = "str_sub(id_hash, 3, 6)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "hash18",
    id_col = "type",
    ids = "HEK_K562_KG1_THP1",
    code = "str_remove(id_hash, '[\\.-].+')"
  )
)