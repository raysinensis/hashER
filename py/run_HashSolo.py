import argparse
import scanpy.external as sce
import pegasus as pg
import anndata
import numpy as np
import pandas as pd

def version():
    print("0.0.1")

"""
Call HashSolo with api to directly output csv 
"""    

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("-c", "--hto", help="path to HTO count matrix, default ./HTO_counts.csv", default="./HTO_counts.csv", type=str)
    parser.add_argument("-o", "--out", help="output file name, default ./HashSolo_results.csv", default="./HashSolo_results.csv.gz", type=str)
    parser.add_argument("-w", "--whitelist", help="barcode whitelist file name", type=str)
    args = parser.parse_args()
    hto_path = args.hto
    res_path = args.out
    whitelist_path = args.whitelist
    hash = anndata.read_csv(hto_path)
    hash.obs = hash.to_df()
    sce.pp.hashsolo(hto, hto.var_names.tolist())
    
    bc = pd.read_csv(whitelist_path, header = None)[0].tolist()
    hash.obs[["negative_hypothesis_probability","singlet_hypothesis_probability","doublet_hypothesis_probability", "Classification"]].loc[bc,].to_csv(res_path")

if __name__ == '__main__': main()
