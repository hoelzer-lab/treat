/*
* Transrate call
*/
workflow TRANSRATE{
  get:
    ASSEMBLIES
    READS
  main:
    TRANSRATE_RUN(ASSEMBLIES, READS)
  emit:
    TRANSRATE.out
}

process TRANSRATE_RUN {
  label "TRANSRATE"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "transrate_out/assemblies.csv"

  input:
  file(assemblies)
  tuple file(leftReads), file(rightReads)


  output:
  file("transrate_out/assemblies.csv")
  
  shell:
  """
  transrate --left ${leftReads} --right ${rightReads} --reference !{params.transcripts} --threads !{params.threads} --output transrate_out --assembly ${assemblies}'
  """ 
}

/* Comments:
*/