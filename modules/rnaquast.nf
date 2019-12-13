/*
* RNAQUAST
*/

workflow RNAQUAST{
  get:
    // TODO: move the params into the workflow call
  main:
    RNAQUAST_SINGLE(params.assemblies.collect(), params.labels, params.reads, params.genome, params.annotation)
    // RNAQUAST_PAIRED(params.assemblies.collect(), params.labels, params.reads, params.genome, params.annotation)
  emit:
    RNAQUAST_SINGLE.out
}

process RNAQUAST_SINGLE {
  label "RNAQUAST"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "short_report.tsv"

  input:
    file(assemblies)
    val(names)
    file(reads)
    file(reference)
    file(annotation)

  output:
    file("short_report.tsv")
  
  shell:
  """
  rnaQUAST.py --single_reads ${reads} \
            --reference ${reference} \
						--gtf ${annotation} \
						--output_dir . \
						--threads !{params.threads} \
						--transcripts ${assemblies} \
						--labels ${names} \
            --no_plots \
            --prokaryote
  """ 
}

process RNAQUAST_PAIRED {
  label "RNAQUAST"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "short_report.tsv"

  input:
    file(assemblies)
    val(names)
    tuple file(leftReads), file(rightReads)
    file(reference)
    file(annotation)

  output:
    file("short_report.tsv")
  
  shell:
    """
    rnaQUAST.py --left_reads ${leftReads} \
                --right_reads ${rightReads} \
                --reference ${reference} \
                --gtf ${annotation} \
                --output_dir . \
                --threads !{params.threads} \
                --transcripts ${assemblies} \
                --labels ${names} \
                --no_plots \
                --prokaryote
    """ 
}

/* Comments:
rnaQUAST is still on python2 in the rnaQUAST.py executable the shebang is still which defaults to #!/usr/bin/env python, but should now default to #!/usr/bin/env python2, because python3 is default on new systems...
*/