# TREAT: TranscRiptome EvaluATion

_Treat your assemblies well!_

A pipeline combining docker containers and nextflow for the evaluation of (_de
novo_) transcriptome assembly results.

````
nextflow main.nf --assemblies test_data/rna-spades.fasta --reads test_data/eco.fastq --threads 2
````