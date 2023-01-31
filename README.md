# hashER
non-optimal scenarios benchmarking of sc(n/multiome)RNA-seq hashing demultiplex tools

# Goal
Some hashing data has just been rushed to the emergency room. It's clearly not doing well -- the data is messy, maybe not even bimodal! So can it be rescued? Which tool(s), often designed and optimized for idealized situations, should be called upon?

# Functionality
1. Iterate over different demultiplexing tools, given the same input files
2. Summarize accuracy and other metrics, given ground truth and cell type info
3. Generate basic QC for hashing data, and suggest appropriate demultiplexing tool

# Non-optimal scenarios and assessment functions
1. Limited sequencing depth - (`get_saturation`)
2. Tag count distribution not bi-modal - (`test_bimodal`)
3. High background contamination from one tag - (`test_background`)
4. Cell type bias for staining and peak count - (`test_typebias`)

# Library Normalization add-on for `Seurat::HTODemux`
We often observe a cell-type bias for hashtag (feature barcoding) quantification and demultiplex results. In brief, if certain cell types attach to the hashtags more readily than others (more surface targets, large size, etc), the peak intensities and thresholds for these cell types would naturally be different, which affects the threshold calling results of HTODemux. Adding a simple library size normalization to the default CLR transformation in our experience often helps improve demultiplex performance. This can be called with `run_HTODemux_Libnorm` from input files, and `HTODemux_Libnorm` for Seurat objects as input.
