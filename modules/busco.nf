/*
* Benchmark universal orthologous groups w/ BUSCO
*/ 
process BUSCO {
  label 'BUSCO'
  publishDir params.outdir, mode:'copy'

  input:
  set assembly_id, file(assembly)
  val busco_dataset
  val threads

  output:
  set assembly_id, file("run_${assembly_id}/short_summary_${assembly_id}.txt")

  shell:
  '''
  run_busco -i !{assembly} -o !{assembly_id} -l /!{busco_dataset}/ -m tran -c !{threads}
  '''
}


