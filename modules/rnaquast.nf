/*
* rnaquast call
*/
process RNAQUAST {
  label 'RNAQUAST'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/short_report.tsv"

  input:
  set val(name), file(assembly)
  set val(read_id), file(reads)
  file(reference)
  file(annotation)

  output:
  set val(name), file("${name}/short_report.tsv")
  
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

/* Comments:
rnaQUAST is still on python2 in the rnaQUAST.py executable the shebang is still which defaults to #!/usr/bin/env python, but should now default to #!/usr/bin/env python2, because python3 is default on new systems...
*/