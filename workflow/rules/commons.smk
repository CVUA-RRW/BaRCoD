import os


shell.executable("bash")


# Helper functions


def get_taxid_fasta_list(wildcards):
    checkpoint_output = checkpoints.export_sequences.get(**wildcards).output['fastadir']
    return expand("fastadump/{fasta}.fa",
                  fasta=glob_wildcards(os.path.join(checkpoint_output, "{fasta}.fa")).fasta)


def generate_db_name(wildcards=None):
    path, dbname = os.path.split(config["blast_db"])
    path, primername = os.path.split(config["primers"])
    return dbname + '_' + primername.split('.')[0]


# General rules


rule prep_taxonomy: 
    output:
        tax="common/taxonomy.json",
    params:
        nodes=config["nodes_dmp"],
        rankedlineage=config["rankedlineage_dmp"],
        taxid=config["parent_node"],
    message:
        "Preparing taxonomy definitions"
    conda:
        "../envs/taxidtools.yaml"
    script:
        "../scripts/filter_taxonomy.py"


rule copy_taxdb:
    output:
        bti = "barcodes/taxdb.bti",
        btd = "barcodes/taxdb.btd",
    params:
        taxdb = config['taxdb'],
    message:
        "Copying taxonomy definitions"
    shell:
        """
        cp {params.taxdb}/taxdb.bti {output.bti}
        cp {params.taxdb}/taxdb.btd {output.btd}
        """
