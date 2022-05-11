# Dereplication, preclustering and pairwise alignement of barcodes


shell.executable("bash")


# Split sequences by taxid


rule get_listing:
    input:
        DB = expand("barcodes/{name}.{ext}", name=generate_db_name(), 
                    ext=["nto", "ntf", "nsq", "not", "nos", "nog", "nin", "nhr", "ndb"]),
        bti = "barcodes/taxdb.bti",
        btd = "barcodes/taxdb.btd",
    output:
        listing = temp("barcodes/db_listing.txt"),
    params:
        dbname = generate_db_name()
    message: 
        "Exporting barcode listing"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        export BLASTDB="barcodes/"
        blastdbcmd -db {params.dbname} -entry all -outfmt '%a\t%T\t%S' > {output.listing}
        """


checkpoint export_sequences:
    input:
        listing = "barcodes/db_listing.txt",
    output:
        fastadir = directory("fastadump/"),
    params:
        blast_DB = f"barcodes/{generate_db_name()}",
        taxdb = config["taxdb"],
    message: 
        "Exporting sequences"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        mkdir {output.fastadir}
        
        export BLASTDB={params.taxdb}
        
        for tax in $(cut -d$'\t' -f2 {input.listing} | sort -u); do
            # retrieving sequences per taxa
            blastdbcmd -db {params.blast_DB} -taxids $tax -outfmt '%f' \
                > {output.fastadir}/$tax.fa
        done
        """


# Filtering


rule dereplicate:
    input:
        fastas = get_taxid_fasta_list,
    output:
        report = temp("derepdump/dereplication.tsv"),
        fasta = temp("fasta/sequences_derep.fa"),
        tmpfa  = temp("derepdump/derep.fa"),
        tmptab = temp("derepdump/derep.txt"),
    message: 
        "Dereplicating sequences"
    conda:
        "../envs/vsearch.yaml"
    shell:
        """
        for file in {input}; do
            
            vsearch --cluster_fast $file \
                    --id 1 \
                    --iddef 1 \
                    --centroids {output.tmpfa} \
                    --uc {output.tmptab} \
                    --quiet
            
            cat {output.tmpfa} >> {output.fasta}
            
            grep -E '^[S|H]' {output.tmptab} \
                 | cut -d$'\t' -f1,9,10 \
                 >> {output.report}
        done
        """


rule filter_seq:
    input:
        fasta = "fasta/sequences_derep.fa",
    output:
        fasta = temp("fasta/sequences.fa"),
    message: 
        "Filtering ambiguous sequences"
    params:
        max_n = config["max_n"],
        min_length = config["min_length"],
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        cutadapt --max-n {params.max_n} --minimum-length {params.min_length} {input.fasta} > {output.fasta}
        """


rule trim_primers:
    input:
        fasta = "fasta/sequences.fa",
    output:
        fasta = temp("fasta/sequences_trim.fa"),
        report = "reports/trimming_report.txt",
        primers_rc = temp("primers_rc.fa"),
    params:
        primers = config["primers"],
    message: 
        "Trimming primers"
    conda: 
        "../envs/cutadapt.yaml"
    shell:
        """
        seqtk seq -r {params.primers} > {output.primers_rc}
        
        cutadapt -g file:{params.primers} \
                 {input.fasta} 2> {output.report} \
            | cutadapt -a file:{output.primers_rc} - \
                       > {output.fasta} 2>> {output.report}
        """


# Alignements


rule pariwise_alignement:
    input:
        fasta = "fasta/sequences_trim.fa" if config["trim_primers"] == True else "fasta/sequences.fa",
    output:
        table = "pairwise_alignment/table.tsv",
        align = "pairwise_alignment/alignments.txt",
    params:
        id = config["min_identity"],
    threads: 
        workflow.cores
    message: 
        "Performing pairwise global alignment"
    conda:
        "../envs/vsearch.yaml"
    shell:
        """
        vsearch --allpairs_global {input.fasta} \
                --iddef 1 \
                --gapext 2I/2E \
                --gapopen 20I/20E \
                --threads {threads} \
                --id {params.id} \
                --userout {output.table} \
                --userfields query+target+id+mism+gaps+alnlen+qlo+qhi+ql+tlo+thi+tl \
                --alnout {output.align}
                
        sed -i '1 i\query\ttarget\tid\tmismatch\tgaps\taln_length\tqstart\tqend\tqlength\ttstart\ttend\ttlength' {output.align}
        """
