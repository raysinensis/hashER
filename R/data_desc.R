hash_data <- setNames(
  data.frame(matrix(ncol = 3, nrow = 0)), 
  c("sample", "ids", "code")
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "A17k",
    ids = "F_M",
    code = "str_sub(pred, 1, 1)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "A45k",
    ids = "F_M",
    code = "str_sub(pred, 1, 1)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "P1mc",
    ids = "3911_9319",
    code = "str_sub(pred, 3, 6)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "P2a",
    ids = "3911_9319",
    code = "str_sub(pred, 3, 6)"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample = "hash18",
    ids = "HEK_K562_KG1_THP1",
    code = ""
  )
)