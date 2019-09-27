params.outdir = 'treat_results'

/*
* Mapping w/ HISAT2 and preprocessing w/ samtools
*/
process HISAT2 {
  publishDir params.outdir, mode:'copy'

  input:
  set assembly_id, file(assembly)
  set read_id, file(reads)
  val threads

  output:
  set assembly_id, file("${assembly_id}.sorted.bam")
  
  shell:
  '''
  hisat2-build !{assembly} !{assembly_id} 
  hisat2 -x !{assembly_id} -U !{reads} -p !{threads} | samtools view -bS | samtools sort -T tmp --threads !{threads} > !{assembly_id}.sorted.bam
  ''' 
}

/*
* Benchmark universal orthologous groups w/ BUSCO
*/ 
process BUSCO {
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


