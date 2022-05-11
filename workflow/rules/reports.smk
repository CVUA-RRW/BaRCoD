# Reporting stats etc


import pandas as pd


shell.executable("bash")


# Reports for barcode recovery


rule recovery_report:
    input:
        sequences = "db_filtering/seq_table.tsv",
        barcodes =  "primer_blaster/barcode_pos.tsv",
    output:
        report = report("reports/recovery_report.tsv", 
                        caption="../report/recovery_report.rst",
                        category="Barcode recovery"),
    message:
        "Collecting barcode recovery stats"
    shell:
        """
        nseq=$(wc -l {input.sequences})
        ntaxid=$(cut -f2 {input.sequences} | sort | uniq -u | wc -l)
        nbarcode=$(wc -l {input.barcodes})
        ntaxidbarcode=$(cut -f2 {input.barcodes} | sort | uniq -u | wc -l)
        echo "DB entries\tTotal taxids\tBarcodes recovered\tTaxids represented" > {output.report}
        echo $nseq\t$ntaxid\t$nbarcode\t$ntaxidbarcode >>  {output.report}
        """


rule primer_hist:
    input:
        primers =  "primer_blaster/primer_blast.tsv",
    output: 
        primer_hist = "reports/primer_sites_distribution.tsv",
    message:
        "Calculating primer binding site distribution"
    shell:
        """
        echo "count\tprimer" > {output.primer_hist}
        tail -n+2 {input.primers} \
            | cut -f2 \
            | cut -d[ -f1 \
            | sort \
            | uniq -c \
            | awk '{{$1=$1}};1' \
            | tr " " "\t" \
            >> {output.primer_hist}
        """


rule barcode_hist:
    input:  
        barcodes =  "primer_blaster/barcode_pos.tsv",
    output:
        barcode_hist = "reports/taxid_barcodes_distribution.tsv",
    message:
        "Calculating barcode distributions per taxid"
    shell:
        """
        echo "count\ttaxid" > {output.barcode_hist}
        tail -n+2 {input.barcodes} \
            | cut -f2 \
            | sort \
            | uniq -c \
            | awk '{{$1=$1}};1' \
            | tr " " "\t" \
            >> {output.barcode_hist}
        """


rule barcode_size_hist:
    input:  
        barcodes =  "primer_blaster/barcode_pos.tsv",
    output:
        size_hist = "reports/barcode_size_distribution.tsv",
    message:
        "Calculating barcode length distributions"
    shell:
        """
        echo "count\ttaxid" > {output.size_hist}
        tail -n+2 {input.barcodes} \
            | cut -f4 \
            | sort \
            | uniq -c \
            | awk '{{$1=$1}};1' \
            | tr " " "\t" \
            >> {output.size_hist}
        """


rule missing_taxids:
    input:
        sequences = "db_filtering/seq_table.tsv",
        barcodes =  "primer_blaster/barcode_pos.tsv",
    output:
        missing = "reports/taxid_wo_barcode.tsv",
        missing_no_deets = temp("reports/taxid_wo_barcode_nosciname.tsv"),
    message:
        "Finding taxids without barcode"
    shell:
        """
        comm -3 \
            <(cut -f2 {input.sequences} | sort -u) \
            <(cut -f2 {input.barcodes} | sort -u) \
            > {output.missing_no_deets}
        
        grep -f {output.missing_no_deets} \
            <(cut -f1,2 {input.sequences}) \
            > {output.missing}
        """


rule recovery_plots:
    input:
        primers = "primer_blaster/primer_blast.tsv",
        barcodes = "primer_blaster/barcode_pos.tsv",
    output:
        primer_hist_plot = report("reports/primer_sites_distribution.pdf",
                             caption="../report/recovery_primerhist.rst",
                             category="Barcode recovery"),
        barcode_hist_plot = report("reports/taxid_barcodes_distribution.pdf",
                              caption="../report/recovery_barcodehist.rst",
                              category="Barcode recovery"),
        size_hist_plot = report("reports/barcode_size_distribution.pdf",
                           caption="../report/recovery_sizehist.rst",
                           category="Barcode recovery"),
    conda:
        "../envs/tidyr.yaml"
    script:
        "../scripts/recovery_plots.R"


# Reports for clustering


rule derep_stats:
    input:
        derep = "derepdump/dereplication.tsv",
        table = "barcodes/db_listing.txt"
    output:
        derep = report("reports/dereplication.tsv",
                       caption="../report/comparison_derep.rst",
                       category="Barcode Dereplication"),
    message: 
        "Collecting dereplication stats"
    shell:
        """
        join --nocheck-order -1 2 -2 1 -t $'\t' \
             <(sort -k2 {input.derep}) \
             <(sort -k1 {input.table}) \
             | sed -e 's/\tH\t/\thit\t/' -e 's/\tS\t/\tcentroid\t/' \
             | sort -k4n \
             > {output}
        
        sed -i '1 i\seqid\ttype\tcentroid\ttaxid\tname' {output.derep}
    """


rule cluster_size:
    input:
        derep = "reports/dereplication.tsv"
    output: 
        sizes = report("reports/cluster_size.tsv",
                       caption="../report/comparison_cluster_sizes.rst",
                       category="Clusters"),
    message: "Collecting cluster sizes"
    run:
        df = pd.read_csv(input['derep'], sep= '\t')
        taxcount =  df.groupby('taxid').agg('count')['seqid'].rename('tax_size')
        centroids = df[df['type'] == 'centroid']
        centroidSize = df.groupby('centroid').nunique()['seqid']+1
        centroidSize = centroidSize[centroidSize.index != '*'].rename("cluster_size")
        dfout = centroids.set_index('seqid').join(centroidSize).fillna(1)
        dfout = dfout.join(taxcount, on = 'taxid').drop(["type", "centroid"], axis = 1)
        dfout = dfout.astype({'cluster_size' : 'int32', 'tax_size': 'int32'})
        dfout["size"] = dfout["cluster_size"].map(str) + "/" + dfout["tax_size"].map(str)
        dfout["rel_cluster_size"] = round(dfout["cluster_size"] / dfout["tax_size"] *100, 2)
        dfout.to_csv(output['sizes'], sep = '\t')


rule get_seq_sizes: 
    input:
        raw = "fasta/sequences.fa",
        trimmed = "fasta/sequences_trim.fa" if config["trim_primers"] == True else "fasta/sequences.fa"
    output:
        raw = temp("reports/fasta_length_raw.tsv"),
        length_trim = temp("reports/fasta_length_trim.txt"),
    params:
        trim = config["trim_primers"]
    message: 
        "Getting barcode length distribution"
    shell:
        """
        # Get seqid length table
        cat {input.raw} \
            | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' \
            > {output.raw}
        
        if [ {params.trim} = "True" ]; then
            # get trimmed lengths
            cat {input.trimmed} \
                | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' \
                > {output.length_trim}
        
        else
            # if not trimming, just merge infos
            touch {output.length_trim} # otherwise snakemake complains about missing output
        
        fi
        """


rule cluster_plot:
    input:
        sizes = "reports/cluster_size.tsv",
    output:
        size_plot = report("reports/cluster_size_distribution.pdf",
                             caption="../report/clustersize_hist.rst",
                             category="Clusters"),
        taxid_size = report("reports/cluster_taxid_distribution.pdf",
                             caption="../report/taxidsize_hist.rst",
                             category="Clusters"),
    message:
        "Plotting cluster size distribution"
    conda:
        "../envs/tidyr.yaml"
    script:
        "../scripts/cluster_plots.R"


# Result tables


rule pairwise_dist:
    input:
        aln = "pairwise_alignment/table.tsv",
        info = "barcodes/db_listing.txt",
        sizes = "reports/cluster_size.tsv"
    output:
        dist = report("reports/distances.tsv",
                      caption="../report/comparison_pairwise_dist.rst",
                      category="Pairwise-distances"),
    message: 
        "Collecting alignement informations"
    run:
        dfaln = pd.read_csv(input.aln, sep='\t', names= ["query", "target", "id", "mismatch", "gaps","aln_length", "qstart", "qend", "qlength", "tstart", "tend", "tlength"])
        dfinfo = pd.read_csv(input.info, sep = '\t', names = ["seqid", "taxid", "name"])
        dfsize = pd.read_csv(input.sizes, sep = '\t')
        dfout = pd.DataFrame({"query": dfaln["query"],
                              "target": dfaln["target"],
                              "distance": dfaln["mismatch"] + dfaln["gaps"]})
        dfout = dfout.join(dfinfo.set_index("seqid"), 
                            on = "query", how = "inner").rename(columns={'taxid' : 'query_taxid', 
                                                                         'name' : 'query_name'})
        dfout = dfout.join(dfinfo.set_index("seqid"), 
                            on = "target", how = "inner").rename(columns={'taxid' : 'target_taxid', 
                                                                         'name' : 'target_name'})
        dfout = dfout.join(dfsize.set_index("seqid")[['size', 'rel_cluster_size']], 
                            on = "query", how = "inner").rename(columns={'size' : 'query_size', 
                                                                         'rel_cluster_size' : 'query_relsize'})
        dfout = dfout.join(dfsize.set_index("seqid")[['size', 'rel_cluster_size']], 
                            on = "target", how = "inner").rename(columns={'size' : 'target_size', 
                                                                         'rel_cluster_size' : 'target_relsize'})
        dfout = dfout[['query', 'query_taxid', 'query_name', 'query_size', 'query_relsize',
                        'target', 'target_taxid', 'target_name', 'target_size', 'target_relsize',
                        'distance']]
        dfout.to_csv(output['dist'], sep='\t', index=False)


rule seq_size_table:
    input:
        raw = "reports/fasta_length_raw.tsv",
        trim = "reports/fasta_length_trim.txt",
        info = "barcodes/db_listing.txt"
    output:
        lengths = report("reports/sequence_lengths.txt",
                         caption="../report/comparison_seqlength.rst",
                         category="Barcode lengths"),
    params:
        trim = config["trim_primers"]
    message: "Formatting sequence length table"
    run:
        dfinfo = pd.read_csv(input.info, sep = '\t', names = ["seqid", "taxid", "name"])
        dfraw = pd.read_csv(input.raw, sep = '\t', names = ["seqid", "length"])
        dfoutraw = dfraw.join(dfinfo.set_index("seqid"),
                              on = "seqid", how = "inner")
        if params.trim:
            dfout = dfoutraw.rename(columns={'length': 'db_length'})
            dftrim = pd.read_csv(input.trim, sep = '\t', names = ["seqid", "length"])
            dfout = dfout.join(dftrim.set_index('seqid'),
                               on = 'seqid', how = 'inner').rename(columns = {'length' : 'trim_length'})
            dfout = dfout[['seqid', 'taxid', 'name', 'db_length', 'trim_length']]
            dfout.to_csv(output[0], sep = '\t', index = False)
        else:
            dfoutraw = dfoutraw[['seqid', 'taxid', 'name', 'length']]
            dfoutraw.to_csv(output['lengths'], sep = '\t', index = False)


rule db_stats:
    input:
        unfiltered = "fasta/sequences_derep.fa",
        seq = "fasta/sequences.fa",
        table = "barcodes/db_listing.txt",
    output:
        table = report("reports/barcode_stats.tsv",
                       caption="../report/comparison_dbinfo.rst",
                       category="Barcodes stats"),
    shell:
        """
        ntaxids=$(cut -d$'\t' -f2 {input.table} | sort -u | wc -l | cut -d$' ' -f1)
        nseq=$(wc -l {input.table} | cut -d$' ' -f1)
        dereps=$(grep -c "^>" {input.seq})
        highN=$(( $(grep -c "^>" {input.unfiltered}) - $(grep -c "^>" {input.seq}) ))
        
        echo "Taxid number\tSequence number\tDereplicated sequences\tLow quality sequences" > {output.table}
        echo $ntaxids\t$nseq\t$dereps\t$highN >> {output.table}
        """


rule get_consensus_level:
    input:
        distance_table = "reports/distances.tsv",
        tax = "common/taxonomy.json",
    output:
        cons = report("reports/consensus.tsv",
                      caption="../report/comparison_consensus.rst",
                      category="Consensus ranks determination"),
    message:
        "Determining consensus ranks"
    params:
        cons_level = config['consensus_level'],
    conda:
        "../envs/taxidtools.yaml"
    script:
        "../scripts/consensus_levels.py"


rule comparison_plots:
    input:
        distance_table = "reports/distances.tsv",
        cons = "reports/consensus.tsv",
    output:
        dist_plot = report("reports/barcodes_distance_plot.pdf",
                      caption="../report/barcodes_distance_plot.rst",
                      category="Pairwise-distances"),
        consensus_plot = report("reports/consensus_plot.pdf",
                      caption="../report/consensus_plot.rst",
                      category="Consensus ranks determination"),
    message:
        "Plotting distances"
    conda:
        "../envs/tidyr.yaml"
    script:
        "../scripts/dist_plots.R"
