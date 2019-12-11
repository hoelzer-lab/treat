workflow Ex90N50{
    main:
        SALMON_INDEX(params.assemblies)
        SALMON_QUASIMAPPING(params.reads, SALMON_INDEX.out)
        ABUNDANCE_ESTIMATES_TO_MATRIX(SALMON_QUASIMAPPING.out)
        CONTIG_EXN50_STATISTIC(params.assemblies, ABUNDANCE_ESTIMATES_TO_MATRIX.out)
        //ALIGN_AND_ESTIMATE_ABUNDANCE(params.assemblies, params.reads)
        //ABUNDANCE_ESTIMATES_TO_MATRIX(ALIGN_AND_ESTIMATE_ABUNDANCE.out)
        //CONTIG_EXN50_STATISTIC()
}

process SALMON_INDEX {
  label 'Ex90N50'

  input:
  tuple val(name), file(assembly)

  output:
  file('*')

  shell:
    """
    salmon index --keepDuplicates -t ${assembly} -i salmon -p !{params.threads}
    """
}

process SALMON_QUASIMAPPING {
  label 'Ex90N50'

  input:
  file(reads)
  file('*')

  output:
  file('*')

  shell:
    """
    salmon quant -l A -i salmon -r ${reads} -o salmon_quant -p !{params.threads}
    """
}

process ABUNDANCE_ESTIMATES_TO_MATRIX {
  label 'Ex90N50'

  input:
  file('*')

  output:
  file('*')

  shell:
  """
  abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --out_prefix salmon_quant --name_sample_by_basedir salmon_quant/quant.sf
  """ 
  }

process CONTIG_EXN50_STATISTIC {
  label 'Ex90N50'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "ExN50.stats"

  input:
  tuple val(name), file(assembly)
  file('*')

  output:
  file('ExN50.stats')
  
  shell:
  """
  contig_ExN50_statistic.pl salmon_quant.isoform.TPM.not_cross_norm ${assembly} > ExN50.stats
  """ 
  }