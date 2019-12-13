/*
* Mapping w/ HISAT2 and preprocessing w/ samtools
*/

workflow MAPPING_SINGLE{
    main:
        HISAT2_SINGLE(params.assemblies, params.reads)  
        SAMTOOLS_FLAGSTATS(HISAT2_SINGLE.out)
        CALCULATE_MAPPED_READS(SAMTOOLS_FLAGSTATS.out)
    emit:
        CALCULATE_MAPPED_READS.out
}

workflow MAPPING_PAIRED{
    main:
        HISAT2_PAIRED(ASSEMBLY)  
        SAMTOOLS_FLAGSTATS(HISAT2_PAIRED.out)
        CALCULATE_MAPPED_READS(SAMTOOLS_FLAGSTATS.out)
    emit:
        CALCULATE_MAPPED_READS.out
}

process HISAT2_SINGLE {
  label 'HISAT2'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}.sorted.bam"

  input:
  tuple val(name), file(assembly)
  file(reads)

  output:
  tuple val(name), file("${name}.sorted.bam")
  
  shell:
  """
  hisat2-build ${assembly} ${name} 
  hisat2 -x ${name} -U ${reads} -p !{params.threads} | samtools view -bS | samtools sort -T tmp --threads !{params.threads} > ${name}.sorted.bam
  """ 
}

process HISAT2_PAIRED {
  label 'HISAT2'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}.sorted.bam"

  input:
  tuple val(name), file(assembly), file(reads)
  //tuple val(read_id), file(reads)

  output:
  set val(name), file("${name}.sorted.bam")
  
  shell:
  """
  hisat2-build ${assembly} ${name} 
  hisat2 -x !{name} -1 ${reads[0]} -2 ${reads[1]} -p !{params.threads} | samtools view -bS | samtools sort -T tmp --threads !{params.threads} > ${name}.sorted.bam
  """ 
}

process SAMTOOLS_FLAGSTATS {
  input:
    tuple val(name), file(assembly)

  output:
    tuple val (name), file('flagstats.txt')
  
  shell:
  """
  samtools flagstat ${assembly} > "flagstats.txt"
  """
}

process CALCULATE_MAPPED_READS {
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}_mapping_stats.txt"

  input:
    tuple val(name), file(mapping_stats)

  output:
    file("${name}_mapping_stats.txt")  
  
  """
  mapping_percentage.py --in ${mapping_stats} --out ${name}_mapping_stats.txt
  """
}

/* Comments:
*/