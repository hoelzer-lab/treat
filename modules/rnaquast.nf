/*
* RNAQUAST
*/

workflow RNAQUAST{
    main:
        RNAQUAST_SINGLE(params.assemblies, params.reads, params.genome, params.annotation)
        // RNAQUAST_SINGLE(assemblies_ch, reads_ch, reference_ch, annotation_ch)

}

process RNAQUAST_SINGLE {
  label 'RNAQUAST'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/short_report.tsv"

  input:
  tuple val(name), file(assembly)
  tuple val(read_id), file(reads)
  file(reference)
  file(annotation)

  output:
  tuple val(name), file("${name}/short_report.tsv")
  
  shell:
  """
  rnaQUAST.py --single_reads ${reads} \
            --reference ${reference} \
						--gtf ${annotation} \
						--output_dir ${name} \
						--threads !{params.threads} \
						--transcripts ${assembly} \
						-l !{name} \
            --prokaryote
  """ 
}

process RNAQUAST_PAIRED {
  label 'RNAQUAST'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/short_report.tsv"

  input:
  tuple val(name), file(assembly)
  tuple val(read_id), file(reads)
  file(reference)
  file(annotation)

  output:
  tuple val(name), file("${name}/short_report.tsv")
  
  shell:
  """
  rnaQUAST.py --left_reads {input.R1} \
						  --right_reads {input.R2} \
						  --reference {input.refGenome} \
						  --gtf {input.refGTF} \
						  -output_dir {params.outDir} \
						  --threads {threads} \
						  --transcripts {params.paths} \
						  --labels {params.names} \
              {params.prokaryote}
  """ 
}




/* Comments:
rnaQUAST is still on python2 in the rnaQUAST.py executable the shebang is still which defaults to #!/usr/bin/env python, but should now default to #!/usr/bin/env python2, because python3 is default on new systems...
*/