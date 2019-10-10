/*
* Mapping w/ HISAT2 and preprocessing w/ samtools
*/
process HISAT2 {
  label 'HISAT2'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}.sorted.bam"

  input:
  set val(name), file(assembly)
  set val(read_id), file(reads)

  output:
  set val(name), file("${name}.sorted.bam")
  
  shell:
  '''
  hisat2-build !{assembly} !{name} 
  hisat2 -x !{name} -U !{reads} -p !{params.threads} | samtools view -bS | samtools sort -T tmp --threads !{params.threads} > !{name}.sorted.bam
  ''' 
}

/* Comments:
*/