run_demuxEM <- function(raw_counts = NULL,
                        raw_h5 = NULL,
                        out_name = "demuxEM.csv.gz",
                        thread = 4) {
  com1 <- paste0(
    "conda run -n base39_0322 ",
    "python ",
    "/nfs/home/rfu/projects/hashER/py/run_demuxEM.py",
    " -d ",
    raw_h5,
    " -c ",
    raw_counts,
    " -o ",
    out_name,
    " -t ",
    thread)
  
  system(com1)
}