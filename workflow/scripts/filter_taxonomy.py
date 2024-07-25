#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import taxidTools as txd


def main(nodes, lineage, merged, taxid, out):
    tax = txd.read_taxdump(nodes, lineage, merged)
    tax.prune(taxid)
    tax.write(out)


if __name__ == '__main__':
    main(snakemake.params['nodes'],
         snakemake.params['rankedlineage'],
         snakemake.params['merged'],
         snakemake.params['taxid'],
         snakemake.output['tax'])