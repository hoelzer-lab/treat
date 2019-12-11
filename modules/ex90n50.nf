workflow EX90N50{
    main:
        SALMON_INDEX(params.assemblies)
        SALMON_QUASIMAPPING(params.reads, SALMON_INDEX.out)
        ABUNDANCE_ESTIMATES_TO_MATRIX(SALMON_QUASIMAPPING.out)
        CONTIG_EXN50_STATISTIC(params.assemblies, ABUNDANCE_ESTIMATES_TO_MATRIX.out)
    emit:
      CONTIG_EXN50_STATISTIC.out
}

process SALMON_INDEX {
  label 'Ex90N50'

  input:
  tuple val(name), file(assembly)

  output:
  file("salmon_${name}/")
  val(name)

  shell:
    """
    salmon index --keepDuplicates -t ${assembly} -i salmon_${name} -p !{params.threads}
    """
}

process SALMON_QUASIMAPPING {
  label 'Ex90N50'

  input:
  file(reads)
  file('*')
  val(name)

  output:
  file('*')
  val(name)

  shell:
    """
    salmon quant -l A -i salmon_${name} -r ${reads} -o salmon_quant_${name} -p !{params.threads}
    """
}

process ABUNDANCE_ESTIMATES_TO_MATRIX {
  label 'Ex90N50'

  input:
  file('*')
  val(name)

  output:
  file('*')
  val(name)

  shell:
  """
  abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --out_prefix salmon_quant_${name} --name_sample_by_basedir salmon_quant_${name}/quant.sf
  """ 
  }

process CONTIG_EXN50_STATISTIC {
  label 'Ex90N50'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}_ExN50.stats"

  input:
  tuple val(name), file(assembly)
  file('*')
  val(name)

  output:
  file("${name}_ExN50.stats")
  
  shell:
  """
  contig_ExN50_statistic.pl salmon_quant_${name}.isoform.TPM.not_cross_norm ${assembly} > ${name}_ExN50.stats
  """ 
  }