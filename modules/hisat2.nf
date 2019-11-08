/*
* Mapping w/ HISAT2 and preprocessing w/ samtools
*/

workflow HISAT2{
    main:
      // if( params.paired )
      //   HISAT2_PAIRED
      // else
        HISAT2_SINGLE(params.assemblies, params.reads)    
}

process HISAT2_SINGLE {
  label 'HISAT2'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}.sorted.bam"

  input:
  tuple val(name), file(assembly)
  tuple val(read_id), file(reads)

  output:
  set val(name), file("${name}.sorted.bam")
  
  shell:
  '''
  hisat2-build !{assembly} !{name} 
  hisat2 -x !{name} -U !{reads} -p !{params.threads} | samtools view -bS | samtools sort -T tmp --threads !{params.threads} > !{name}.sorted.bam
  ''' 
}

process HISAT2_PAIRED {
  label 'HISAT2'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}.sorted.bam"

  input:
  tuple val(name), file(assembly)
  tuple val(read_id), file(reads)

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