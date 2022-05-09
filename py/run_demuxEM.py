import argparse
import demuxEM
import pegasus as pg
import anndata
import pegasusio as io
from itertools import compress

def version():
    print("0.0.1")

"""
Call demuxEM with api to directly output csv 
"""    

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-d", "--h5", help="path to raw feature matrix H5 file from cellranger, default ./raw_feature_bc_matrix.h5", default="./raw_feature_bc_matrix.h5", type=str)
    parser.add_argument("-c", "--hto", help="path to HTO count matrix, default ./HTO_counts.csv", default="./HTO_counts.csv", type=str)
    parser.add_argument("-o", "--out", help="output file name, default ./demuxDE_results.csv", default="./demuxDE_results.csv.gz", type=str)
    parser.add_argument("-t", "--thread", help="number of threads", default=1, type=int)
    parser.add_argument("-n", "--min", help="minimum UMIs to keep barcode", default=100, type=int)
    args=parser.parse_args()
    drop_path = args.h5
    hto_path = args.hto
    res_path = args.out
    nthread = args.thread
    nmin = args.min
    dat = io.read_input(drop_path)
    hash = io.MultimodalData(anndata.read_csv(hto_path))
    demuxEM.estimate_background_probs(hash)
    demuxEM.demultiplex(dat, hash, min_signal=10.0, alpha=0.0, alpha_noise=1.0, tol=1e-06, n_threads=nthread)
    
    # out, filter to 100UMI drops
    umis = (dat.X.sum(axis = 1) >= nmin).tolist()
    umis_ind = [item for sublist in umis for item in sublist]
    dat.obs[["demux_type", "assignment"]].loc[list(compress(dat.obs_names.tolist(), umis_ind)), ].to_csv(res_path)
if __name__ == '__main__': main()
