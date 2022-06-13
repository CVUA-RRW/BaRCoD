# BaRCoD - Barcode Recovery and Comparison from Database

A snakemake workflow aiming at recovering and analyzing barcodes for metabarcoding experiments.
Using a set of primers, finds possible amplicon in the database (or a taxonomic subset thereof) and performs pairwise alignements of barcodes.

Orignally RRW-Primerblast and BAnalyzer, now united.

## Usage

Download the repository and create a conda environment containing snakemake.

```bash
conda create -c bioconda -c conda-forge --name snakemake snakemake
```

Then modify the config file template in `config/config.yaml` to suit your needs
and run the workflow with:

```bash
snakemake --use-conda --cores 1 --configfile /path/to/config.yaml
```

More options are available, see the [snakemake's documentation](https://snakemake.readthedocs.io/en/stable/).

## Sequence database

To run the pipeline you will need to provide a BLAST-formated reference sequence database.
If you already have a fasta file with your sequences follow the [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279688/)
to know how to format it.

You will also need to provide the [taxdb](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz) and 
[taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)
files available from the NCBI server.

## Workflow

### Barcode recovery

First the provided primers will be used to scan the (BLAST-formatted) database and find putative 
primer binding sites. Possible amplicons will then be extracted between consecutive primer matches:

```
-> indicates a match on the '+' strand 
<- indicates a match on the '-' strand

primer matches:      <- ->  ->    <-  ->
database sequence: =========================
flanked sequence:           ========

primer matches:      -> <-  ->    <-  ->
database sequence: =========================
flanked sequence:    =====  ========

primer matches:      <- <-  ->    ->  ->
database sequence: =========================
flanked sequence:  
```

### Pairwise comparison

Sequences within a Taxonomic node will be clustered prior to the alignement and calculation of the Hamming distance.
This allows to reduce redundancy in the database. Dereplication is performed by grouping 
sequences with 100% identity **within a taxonomic node**. 
Note that alignement of ambiguous nucleotides **never** incures a penalty. Therefore sequences
containing ambiguous nuclotides can be clustered with sequences that contain strict nucleotide (A, T, U, C, G). 

Then each sequence cluster will be compaired in pairs to determine the pair-wise Hamming distance.
Note that only pairs with at least the identity level defined in the config file be compaired.

If no distance is available for a given taxid it may be that:
* there is no barcode for this taxid
* there is no other barcode within the identtity threshold

To differentiate between these two option check the output of the Barcode recovery process.

### Consensus determination

For each cluster a consensus rank will be determined. This reflects the expected consensus
that can be reached for this sequence considering all sequences within the given identtity 
threshold and with a given agreement level.

## Credits

This would not be possible without: 

* Snakemake
* BLAST+
* Biopython
* Cutadapt
* TaxidTools
* VSEARCH
* R tidyverse

## Contributing

For new features or to report bugs please submit issues directly on the online repository.

## License

This project is licensed under a BSD 3-Clauses License, see the LICENSE file for details.

## Author

For questions about the pipeline, problems, suggestions or requests, feel free to contact:

Grégoire Denay, Chemisches- und Veterinär-Untersuchungsamt Rhein-Ruhr-Wupper 

<gregoire.denay@cvua-rrw.de>