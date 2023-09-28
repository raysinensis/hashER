call_XY <- function(so, 
                    genelist1 = c("Ddx3y", "Eif2s3y", "Kdm5d","Uty"), 
                    genelist2 = c("Xist", "Tsix", "Frmd7"),
                    name1 = "M",
                    name2 = "F",
                    namedoublet = "D",
                    namenegative = "N") {
  sotemp <- so
  sotemp <- AddModuleScore(sotemp, list(genelist1), name = "M")
  sotemp <- AddModuleScore(sotemp, list(genelist2), name = "F")
  sotemp$call <- sotemp@meta.data %>% mutate(call = case_when(
    F1 > 0 & M1 < 0 ~ name2,
    F1 < 0 & M1 > 0 ~ name1,
    F1 > 0 & M1 > 0 ~ namedoublet,
    TRUE ~ namenegative
  )) %>% pull(call)
}