# Workflow configuration

To configure the workflow, modifiy `config/config.yaml` according to your needs.

| Parameter | Description |
| --- | --- |
| workdir | PAth to the working (output) directory |
| threads | Number of threds assigned to the workflow |
| blast_db | Path to the BLAST formatted database, including the basename but no extension (e.g. /path/to/BLAST/nt) |
| taxdb | Path to the folder containing the taxdb files (taxdb.bti and taxdb.btd) |
| rankedlineage_dmp | Path to the rankedlineage.dmp file |
| nodes_dmp | Path to the node.dmp file |
| parent_node | A taxid to limit the BLAST search to a given branch (e.g. 7742 to limit the search to Vertebrates) |
| primers | Path to the fasta file containing the primer sequences |
| primerBlast_coverage | Minimal primer sequence coverage (in percent) |
| primerBlast_identity | Minimal identity to the primer query (in percent) |
| trim_primers | Should primers be trimmed before calculating the hamming distances (True or False) |
| max_n | Maximum number of N nulceotide allowed per barcode |
| min_length | Minimal barcode length |
| max_length | Maximal barcode length |
| min_identity | Minimal identity level to calculate hamming distance |
| consensus_level | Proportion of agreement needed to reach a taxonomic consensus (0.51 to 1.0) |