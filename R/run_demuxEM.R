run_demuxEM <- function(raw_counts = NULL,
                        whitelist = NULL,
                        raw_h5 = NULL,
                        out_name = "demuxEM.csv.gz",
                        whitelist_rm = "-1",
                        thread = 4) {
  com1 <- paste0(
    "conda run -n base39_0322 ",
    "python ",
    "/nfs/home/rfu/projects/hashER/py/run_demuxEM.py",
    " -c ",
    raw_counts,
    " -o ",
    out_name,
    " -t ",
    thread)
  
  if (!is.null(raw_h5) & file.exists(raw_h5)) {
    com1 <- paste0(
      com1,
      " -d ",
      raw_h5
    )
  }
  
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