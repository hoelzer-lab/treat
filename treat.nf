#!/usr/bin/env nextflow

/*
* TREAT WORKFLOW -- TanscRiptome EvaluATion
*
* Author: hoelzer.martin@gmail.com
* Author: lasse.faber@gmail.com
*/

nextflow.preview.dsl=2

// terminal prints
println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.threads"
println "Output dir name: $params.output\u001B[0m"
println " "}

if (params.help) { exit 0, helpMSG() }
if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.assemblies == '') {exit 1, "--assemblies is a required parameter"}

// Queue channels (can only be used once)
assemblies_ch = Channel
              .fromPath( params.assemblies.tokenize(',') )
              .flatMap{ files(it) }
              .map{ file -> tuple(file.simpleName, file) }

assemblies_simplename = Channel
              .fromPath(params.assemblies.tokenize(','))
              .map { file -> file.simpleName }

assemblies_rnaquast = Channel
              .fromPath( params.assemblies.tokenize(',') )
              .flatMap{ files(it) }
              .view()
              
assemblies_rnaquast_labels = Channel
              .fromPath(params.assemblies.tokenize(','))
              .map { file -> file.simpleName }
              .collect()
              .map{it -> it.join(' ')}

// Need to be a value channel, because we want to use this unlimited times.
reads_ch = Channel.value( file(params.reads) )

reference_ch = Channel.value( file(params.reference) )

transcripts_ch = Channel.value( file(params.transcripts) )

annotation_ch = Channel.value( file(params.annotation) )

// illumina reads input & --list support. MIGHT be nicer as initial input read in

// if (params.csv) { illumina_input_ch = Channel
//     .fromPath( params.csv, checkIfExists: true )
//     .splitCsv()
//     .map { row -> ["${row[0]}", "${row[1]}", "${row[2]}"] }
//     .branch{
//       paired: it[2] != ""
//       single: it[2] == ""
//     }
//     .set{result}
// }


// illumina_input_ch.paired.subscribe{println "$it"}
// result.paired.view{"$it"}
// result.single.view{"$it"}


// MAIN WORKFLOW

include './modules/hisat2' params(output: params.output, dir: params.mappingdir, threads: params.threads, assemblies: assemblies_ch, reads: reads_ch)
include './modules/busco' params(output: params.output, dir: params.buscodir, threads: params.threads, assemblies: assemblies_ch, busco: params.busco)
include './modules/transrate' params(output: params.output, dir: params.transratedir, threads: params.threads)
include './modules/rnaquast' params(output: params.output, dir: params.rnaquastdir, threads: params.threads, genome: reference_ch, annotation: annotation_ch, assemblies: assemblies_rnaquast, reads: reads_ch, labels: assemblies_rnaquast_labels)
include './modules/detonate' params(output: params.output, dir: params.detonatedir, threads: params.threads, reference: params.reference, transcripts: transcripts_ch, reads: reads_ch, assemblies: assemblies_ch)
include './modules/ex90n50' params(output: params.output, dir: params.ex90n50dir, threads: params.threads, reads: reads_ch, assemblies: assemblies_ch)
include './modules/extract_metrics' params(assemblies: assemblies_simplename, output: params.output, detonate_dir: params.detonatedir)


// TRANSRATE(assemblies_ch)

workflow {
    main:
        HISAT2()
        BUSCO()
        RNAQUAST()
        DETONATE()
        EX90N50()
        HEATMAP(HISAT2.out.collect(), 
            BUSCO.out.collect(), 
            RNAQUAST.out, 
            DETONATE.out.kc.collect(), 
            DETONATE.out.contig.collect(), 
            DETONATE.out.rsem.collect(), 
            EX90N50.out.collect())
}


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
    nextflow run hoelzer/treat --assemblies test_data/rna-spades.fasta --reads test_data/eco.fastq --threads 4 --busco firmicutes_odb9

    ${c_yellow}Mandatory:${c_reset}
    ${c_green}--assemblies${c_reset}    e.g.: 'trinity.fasta spades.fasta orp.fasta' or '*.fasta' or '*/*.fasta'
    ${c_green}--reads${c_reset}         e.g.: 'trinity.fastq' or '*.fastq' or '*/*.fastq'
    ${c_green}--reference${c_reset}     reference genome
    ${c_green}--annotation${c_reset}    annotation file in gtf format corresponding to the reference file
    ${c_green}--busco${c_reset}         the database used with BUSCO, see https://busco.ezlab.org/v2/frame_wget.html for a full list of available data sets and select one [default: $params.busco]

    ${c_green}--csv${c_reset}     csv input


    ${c_yellow}Options${c_reset}
    --threads                max cores for local use [default: $params.threads]
    --mem                    max memory in GB for local use [default: $params.mem]
    --output                 name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 standard [uses docker], conda, googlegenomics [default: standard] ${c_reset}
    """.stripIndent()
}
