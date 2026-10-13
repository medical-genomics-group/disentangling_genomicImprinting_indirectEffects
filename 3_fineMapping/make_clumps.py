#%%
import pandas as pd
import numpy as np
import networkx as nx
import argparse
import csv

def main(r2,ld_file,path_out):
    ld=pd.read_csv(ld_file, sep=r"\s+",usecols=["SNP_A","SNP_B","R2"])
    ld=ld[ld["R2"]>=r2]

    G = nx.Graph()
    for _, row in ld.iterrows():
        G.add_edge(row["SNP_A"], row["SNP_B"])

    clusters=list(nx.connected_components(G))
    with open(path_out, "w") as f:
        wr = csv.writer(f)
        wr.writerows(clusters)
    print("Saved file")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Clump in groups from r2 plink .ld file')
    parser.add_argument('--r2', type=float, default=0.5, help='r2 threshold to group SNPs, default=0.5')
    parser.add_argument('--ld_file', type=str, help='.ld file with r2 from plink -r2', required = True)
    parser.add_argument('--path_out', type=str, required = True)

    args = parser.parse_args()
    print(args)

    main(r2=args.r2,ld_file=args.ld_file,path_out=args.path_out)
    print("Done")

