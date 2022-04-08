#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import taxidTools as txd


def get_consensus_rank(tax, idlist, cons_level):
    return tax.consensus([str(t) for t in idlist], cons_level).rank


def main(table, taxonomy, outfile, cons_level):
    df = pd.read_csv(table, sep = '\t')
    
    df["entries"] = df["query"] + df["query_size"]
    
    dfout = pd.DataFrame()
    tax = txd.load(taxonomy)
    
    for entry in set(df["entries"]):
        sub = df[df["entries"] == entry].reindex()
        
        new = sub.head(1)[["query_name", "query_taxid", "query_size", "query_relsize"]]
        
        for i in range(0,4):
            taxids = list(sub[sub["distance"] <= i]["target_taxid"])
            taxids.extend(list(sub.head(1)["query_taxid"]))
            rank = get_consensus_rank(tax, taxids, cons_level)
            
            new[f"Consensus rank with {i} mismatches"] = rank
        
        dfout = pd.concat([dfout, new])
        
        dfout.to_csv(outfile, sep = '\t')


if __name__ == '__main__':
    main(snakemake.input['distance_table'],
         snakemake.input['tax'],
         snakemake.output['cons'],
         snakemake.params['cons_level'])