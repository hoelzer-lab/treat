#!/usr/bin/env nextflow

/*
* TREAT WORKFLOW -- TanscRiptome EvaluATion
*
* Author: hoelzer.martin@gmail.com
*/

nextflow.preview.dsl=2
include 'modules/treat'

def helpMSG() {
    log.info """

    Usage:
    nextflow run hoelzer/treat --assemblies test_data/rna-spades.fasta --reads test_data/eco.fastq --threads 4 

    Mandatory:
    --assemblies    e.g.: 'trinity.fasta spades.fasta orp.fasta' or '*.fasta' or '*/*.fasta'
    --reads         e.g.: 'trinity.fastq' or '*.fastq' or '*/*.fastq'
    --reference     reference genome
    --annotation    annotation file in gtf format corresponding to the reference file

    Options
    --threads                max cores for local use [default $params.threads]
    --mem                    max memory in GB for local use [default $params.mem]
    --output                 name of the result folder [default $params.output]

    -with-report rep.html    gives a detailed CPU and RAM usage report in rep.html

    """.stripIndent()
}

if (params.help) { exit 0, helpMSG() }
if (params.assemblies == '') {exit 1, "--assemblies is a required parameter"}

// file channels
assemblies_ch = Channel
              .fromPath(params.assemblies)
              .map { file -> tuple(file.simpleName, file) }

reads_ch = Channel
              .fromPath(params.reads)
              .map { file -> tuple(file.simpleName, file) }


HISAT2( assemblies_ch, reads_ch, params.threads)
BUSCO( assemblies_ch, params.busco_dataset, params.threads )