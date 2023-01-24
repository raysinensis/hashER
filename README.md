# hashER
non-optimal scenarios benchmarking of sc(n/multiome)RNA-seq hashing demultiplex tools

# Goal
Your new hashing data has just been rushed to the emergency room. Visually, it's clearly not doing well -- the data is messy, maybe not even bimodal! So can it be rescued? Which tool(s), often designed and optimized for idealized situations, should be called upon?

# Functionality
1. Iterate over different demultiplexing tools, given the same input files
2. Summarize accuracy and other metrics, given ground truth and cell type info
3. Generate basic QC for hashing data, and suggest appropriate demultiplexing tool
