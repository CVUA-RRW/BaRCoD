# Generates a BLAST formatted database of predicted amplicons from an existing
# BLAST database and a fasta file containing primers.


shell.executable("bash")


# Filtering Input database


rule get_taxid_from_db:
    output:
        taxlist = temp("db_filtering/taxid_list.txt"),
    params:
        blast_db = config["blast_db"],
        taxdb = config["taxdb"],
    message: 
        "Extracting all Taxid from database"
    conda: 
        "../envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastdbcmd -db {params.blast_db} -tax_info -outfmt '%T' > {output.taxlist}
        """


rule filter_taxid:
    input:
        taxlist = "db_filtering/taxid_list.txt",
        tax = "common/taxonomy.json",
    output:
        mask = "db_filtering/taxid_mask.txt",
    params:
        taxid = config["parent_node"],
    message: 
        "Limiting the taxonomy to a single branch"
    conda:
        "../envs/taxidtools.yaml"
    script:
        "../scripts/make_blast_mask.py"


rule get_seqidlist:
    input:
        mask = "db_filtering/taxid_mask.txt",
    output:
        seqids = temp("db_filtering/seqids.txt"),
        binary = temp("db_filtering/seqids.acc"),
        id_table = temp("db_filtering/seq_table.tsv"),
    params:
        blast_db  =config["blast_db"],
        taxdb = config["taxdb"],
    message: 
        "Listing sequences to screen"
    conda: 
        "../envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastdbcmd -db {params.blast_db} -taxidlist {input.mask} -outfmt '%a\t%T\t%S' > {output.id_table}
        
        cat {output.id_table} | cut -d$'\t' -f1 > {output.seqids}
        
        blastdb_aliastool -seqid_file_in {output.seqids} -seqid_file_out {output.binary}
        """


# Primer BLASTer


rule primers_explicit:
    output:
        primers = "primer_blaster/primers_explicit.fa",
    params:
        primers = config["primers"],
    message:
        "Making primer sequences explicit"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/IUPAC_translate.py"


rule find_primer_matches:
    input:
        seqids = "db_filtering/seqids.txt",
        binary = "db_filtering/seqids.acc",
        primers = "primer_blaster/primers_explicit.fa",
    output:
        blast = "primer_blaster/primer_blast.tsv",
    params:
        blast_db = config["blast_db"],
        taxdb = config["taxdb"],
        cov = config["primerBlast_coverage"],
        identity = config["primerBlast_identity"],
    threads: 
        config["threads"]
    message:
        "Primer BLASter"
    conda: 
        "../envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        blastn -db {params.blast_db} \
            -query {input.primers} \
            -task blastn-short \
            -seqidlist {input.binary} \
            -outfmt '6 saccver qseqid staxid sstart send length sstrand mismatch' \
            -ungapped -qcov_hsp_perc {params.cov} -perc_identity {params.identity} \
            -subject_besthit \
            -max_target_seqs  1000000000 \
            -num_threads {threads} \
                | sort -k1 | sed '1 i\seqid\tquery\ttaxid\tstart\tend\tlength\tstrand\tmismatch' > {output.blast}
        """


rule extract_barcodes:
    input:
        blast = "primer_blaster/primer_blast.tsv",
    output:
        positions = "primer_blaster/barcode_pos.tsv",
    message: 
        "Extracting barcodes coordinates"
    script:
        "../scripts/extract_barcodes.py"


rule extract_barcodes_seq:
    input:
        positions = "primer_blaster/barcode_pos.tsv",
    output:
        sequences = "barcodes/{name}.fasta",
    message: 
        "Extracting barcode sequences"
    params:
        blast_db = config["blast_db"],
        taxdb = config["taxdb"],
    conda: 
        "../envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        while IFS=$'\t' read -r acc tax start stop length; do
            blastdbcmd -entry $acc \
                -db {params.blast_db} \
                -range $start-$stop \
                -outfmt %f | sed -e 's/:[[:digit:]-]*//' >> {output.sequences} &
        done < {input.positions}
        """


# Making a BLAST database


rule make_barcode_db:
    input: 
        fasta = f"barcodes/{generate_db_name()}.fasta",
        id_table = "db_filtering/seq_table.tsv",
    output:
        taxid_mapper = temp("barcodes/taxmap.tsv"),
        DB = expand("barcodes/{name}.{ext}", name=generate_db_name(), 
                    ext=["nto", "ntf", "nsq", "not", "nos", "nog", "nin", "nhr", "ndb"]),
    params: 
        blast_db = config["blast_db"],
        taxdb = config["taxdb"],
        dbname = f"barcodes/{generate_db_name()}",
    message: "Formatting barcodes to BLAST database"
    conda: "../envs/blast.yaml"
    shell:
        """
        export BLASTDB={params.taxdb}
        
        cat {input.id_table} | cut -d$'\t' -f1,2 > {output.taxid_mapper}
        
        makeblastdb -in {input.fasta} -dbtype nucl -parse_seqids -blastdb_version 5 -taxid_map {output.taxid_mapper} -out {params.dbname}
        """

