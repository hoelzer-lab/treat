#!/usr/bin/env nextflow

/*
* TREAT WORKFLOW -- TanscRiptome EvaluATion
*
* Author: hoelzer.martin@gmail.com
* Author: lasse.faber@gmail.com
*/

nextflow.preview.dsl=2

if (params.help) { exit 0, helpMSG() }
if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.assemblies == '') {exit 1, "--assemblies is a required parameter"}

// file channels
assemblies_ch = Channel
              .fromPath(params.assemblies)
              .map { file -> tuple(file.simpleName, file) }

reads_ch = Channel
              .fromPath(params.reads)
              .map { file -> tuple(file.simpleName, file) }

// DATABASES

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use.
It also comes with a "auto-download" if a database is not available. Doing it the following way:
1. take userinput DB 2. add a cloud preload DB if cloud profile 3. check if the preload file exists (local or cloud)
4. if nothing is true -> download the DB and store it in the "preload" section (either cloud or local for step 3.)
*/

// get BUSCO db
include 'modules/buscoGetDB' params(db: params.busco_dataset)
buscoGetDB() 
db_busco = buscoGetDB.out


include 'modules/hisat2' params(output: params.output, dir: params.mappingdir, threads: params.threads)
include 'modules/busco' params(output: params.output, db: db_busco, dir: params.buscodir, threads: params.threads)

HISAT2(assemblies_ch, reads_ch)
//BUSCO(assemblies_ch)


def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    TREAT
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run hoelzer/treat --assemblies test_data/rna-spades.fasta --reads test_data/eco.fastq --threads 4 

    ${c_yellow}Mandatory:${c_reset}
    ${c_green}--assemblies${c_reset}    e.g.: 'trinity.fasta spades.fasta orp.fasta' or '*.fasta' or '*/*.fasta'
    ${c_green}--reads${c_reset}         e.g.: 'trinity.fastq' or '*.fastq' or '*/*.fastq'
    ${c_green}--reference${c_reset}     reference genome
    ${c_green}--annotation${c_reset}    annotation file in gtf format corresponding to the reference file

    ${c_yellow}Options${c_reset}
    --threads                max cores for local use [default $params.threads]
    --mem                    max memory in GB for local use [default $params.mem]
    --output                 name of the result folder [default $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 standard, googlegenomics [default: standard] ${c_reset}
    """.stripIndent()
}