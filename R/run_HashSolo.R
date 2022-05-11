run_HashSolo <- function(raw_counts = NULL,
                        whitelist = NULL,
                        out_name = "HashSolo.csv.gz",
                        whitelist_rm = "-1") {
  com1 <- paste0(
    "conda run -n base39_0322 ",
    "python ",
    "/nfs/home/rfu/projects/hashER/py/run_HashSolo.py",
    " -c ",
    raw_counts,
    " -o ",
    out_name)
  
  if (!is.null(whitelist)) {
    com1 <- paste0(
      com1,
      " -w ",
      whitelist
    )
    if (!is.null(whitelist_rm)) {
      com1 <- paste0(
        com1,
        " -r ",
        whitelist_rm
      )
    }
  }
  
  system(com1)
}