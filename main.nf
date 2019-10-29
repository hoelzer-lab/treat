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


assemblies_fullname = Channel
              .fromPath(params.assemblies)
              .map { file -> tuple(file.fileName, file) }

// Need to be a value channel, because we want to use this unlimited times.
reads_ch = Channel.value( tuple( file(params.reads).simpleName, params.reads) )

reads_detonate = Channel.value( file(params.reads) )

reference_ch = Channel.value( file(params.reference) )

transcripts_ch = Channel.value( file(params.transcripts) )

annotation_ch = Channel.value( file(params.annotation) )

// illumina reads input & --list support. MIGHT be nicer as initial input read in
/*
if (params.illumina && params.list) { illumina_input_ch = Channel
    .fromPath( params.illumina, checkIfExists: true )
    .splitCsv()
    .map { row -> ["${row[0]}", [file("${row[1]}"), file("${row[2]}")]] }
    .view() }
else if (params.illumina) { illumina_input_ch = Channel
    .fromFilePairs( params.illumina , checkIfExists: true )
    .view() }
*/

// DATABASES

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use.
It also comes with a "auto-download" if a database is not available. Doing it the following way:
1. take userinput DB 2. add a cloud preload DB if cloud profile 3. check if the preload file exists (local or cloud)
4. if nothing is true -> download the DB and store it in the "preload" section (either cloud or local for step 3.)
*/

// get BUSCO db
// include 'modules/buscoGetDB' //params(db: params.busco)
// buscoGetDB(params.busco) 
// db_busco = buscoGetDB.out


// MAIN WORKFLOW

include './modules/hisat2' params(output: params.output, dir: params.mappingdir, threads: params.threads)
include './modules/busco' params(output: params.output, dir: params.buscodir, threads: params.threads)
include './modules/transrate' params(output: params.output, dir: params.transratedir, threads: params.threads)
include './modules/rnaquast' params(output: params.output, dir: params.rnaquastdir, threads: params.threads, reference: params.reference, annotation: params.annotation)
include './modules/detonate' params(output: params.output, dir: params.detonatedir, threads: params.threads, reference: params.reference, transcripts: params.transcripts)


// HISAT2_SINGLE(assemblies_ch, reads_ch)
// BUSCO(assemblies_ch, db_busco)
// TRANSRATE(assemblies_ch)
// RNAQUAST_SINGLE(assemblies_ch, reads_ch, reference_ch, annotation_ch)

// DETONATE(assemblies_ch, transcripts_ch)
workflow {
    main:
        DETONATE()
}
workflow DETONATE {

  main:
    // performed once   
    ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS(transcripts_ch)
    RSEM_PREPARE_REFERENCE(transcripts_ch)
    RSEM_CALCULATE_EXPRESSION(RSEM_PREPARE_REFERENCE.out, reads_detonate)
    SAMTOOLS_SORT(RSEM_CALCULATE_EXPRESSION.out.filter{ "~/*.bam/" })
    REF_EVAL_ESTIMATE_TRUE_ASSEMBLY(RSEM_PREPARE_REFERENCE.out, RSEM_CALCULATE_EXPRESSION.out, SAMTOOLS_SORT.out)
    RSEM_PREPARE_REFERENCE_2(REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    RSEM_CALCULATE_EXPRESSION_2(RSEM_PREPARE_REFERENCE_2.out, reads_detonate)
    DETONATE7(reads_detonate, SAMTOOLS_SORT.out)
    // performed for each assembly
    RSEM_EVAL_CALCULATE_SCORE(ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS.out, reads_detonate, assemblies_ch)
    REF_EVAL_KC(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out, RSEM_CALCULATE_EXPRESSION_2.out)
    BLAT_A_TO_B(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    BLAT_B_TO_A(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    REF_EVAL_CONTIG(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out, BLAT_A_TO_B.out, BLAT_B_TO_A.out)
  emit:
    RSEM_CALCULATE_EXPRESSION.out
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
