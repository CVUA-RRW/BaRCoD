#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import taxidTools as txd


def get_consensus_rank(tax, idlist, cons_level):
    return tax.consensus([str(t) for t in idlist], cons_level, ignore_missing=True).rank


def main(table, taxonomy, clusters, outfile, cons_level):
    df = pd.read_csv(table, sep = '\t')

    df["entries"] = df["query"] + df["query_size"]

    dfout = pd.DataFrame()
    tax = txd.read_json(taxonomy)

    max_dis = max(df['distance'])+1

    for entry in set(df["entries"]):
        sub = df[df["entries"] == entry].reindex()
        
        new = sub.head(1)[["query_name", "query_taxid", "query_size", "query_relsize"]]
        
        for i in range(0, max_dis):
            taxids = list(sub[sub["distance"] <= i]["target_taxid"])
            taxids.extend(list(sub.head(1)["query_taxid"]))
            rank = get_consensus_rank(tax, taxids, cons_level)
            
            new[f"{i} mismatches"] = rank
        
        dfout = pd.concat([dfout, new])

    # For missing taxids, check if missing or if no match with id threshold
    clus = pd.read_csv(clusters, sep = '\t')
    clus_tax = set(clus['taxid'])
    dfout_tax = set(dfout['query_taxid'])

    missing = clus_tax.difference(dfout_tax)

    clus = clus[clus['taxid'].isin(missing)]

    clus["entries"] = clus["taxid"].astype(str) + clus["size"]

    for entry in set(clus["entries"]):
        sub = clus[clus["entries"] == entry].reindex()
        
        new = sub.head(1)[["name", "taxid", "size", "rel_cluster_size"]]
        new.columns = ["query_name", "query_taxid", "query_size", "query_relsize"]
        
        for i in range(0, max_dis):
            new[f"{i} mismatches"] = tax.getRank(new.iloc[0]['query_taxid'])
        
        dfout = pd.concat([dfout, new])

    dfout.to_csv(outfile, sep = '\t')


if __name__ == '__main__':
    main(snakemake.input['distance_table'],
         snakemake.input['tax'],
         snakemake.input['clusters'],
         snakemake.output['cons'],
         snakemake.params['cons_level'])