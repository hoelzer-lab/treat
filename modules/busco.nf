/*
* Benchmark universal orthologous groups w/ BUSCO
*/ 
process BUSCO {
  label 'BUSCO'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "short_summary_${name}.txt"

  input:
  set val(name), file(assembly)

  output:
  set val(name), file("short_summary_${name}.txt")

  shell:
  '''
  run_busco -i !{assembly} -o !{name} -l !{params.db}/ -m tran -c !{params.threads}
  cp run_!{name}/short_summary_!{name}.txt .
  '''
}

/* Comments:
*/
