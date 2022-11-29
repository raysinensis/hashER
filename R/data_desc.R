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

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "DR",
    id_col = "seurat_clusters",
    ids = "tag1_tag2_tag3_tag4_tag5_tag6_tag7_tag8",
    code = ""
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "dual",
    id_col = "seurat_clusters",
    ids = "BBO21_BBO22_BBO23_BBO25_BBO27_BBO28_BBO29_BBO30_BBO44_BBO51_BBO52_BBO53",
    code = ""
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "amc",
    id_col = "region",
    ids = "MC_A",
    code = "str_remove(id_hash, '_.+')"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "asap",
    id_col = "region",
    ids = "MC_A",
    code = "str_remove(id_hash, '_.+')"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "ag",
    id_col = "orig.ident",
    ids = "MC_A",
    code = "str_remove(id_hash, '_.+')"
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "Landau",
    id_col = "orig.ident",
    ids = "Landau-B1-P2-1_Landau-B1-P2-2_Landau-B1-P2-3_Landau-B1-P2-4",
    code = ""
  )
)

hash_data <- hash_data %>% rbind(
  data.frame(
    sample_id = "SAU",
    id_col = "freemuxlet_id",
    ids = "pt025_pt507_pt620_pt695_pt940",
    code = "str_sub(id_hash, 1, 5)"
  )
)