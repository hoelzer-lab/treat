workflow Ex90N50{
    main:
        ALIGN_AND_ESTIMATE_ABUNDANCE(params.assemblies, params.reads)
        ABUNDANCE_ESTIMATES_TO_MATRIX(ALIGN_AND_ESTIMATE_ABUNDANCE.out)
        //CONTIG_EXN50_STATISTIC()
}

process ALIGN_AND_ESTIMATE_ABUNDANCE {
  label 'Ex90N50'

  input:
  tuple val(name), file(assembly)
  file(reads)

  output:
  file('*')
  
  shell:
  """
  align_and_estimate_abundance.pl --thread_count !{params.threads} --transcripts ${assembly} --seqType fq --single ${reads} --est_method salmon --trinity_mode --prep_reference --output_dir salmon_${name}
  """ 
  }


process ABUNDANCE_ESTIMATES_TO_MATRIX {
  label 'Ex90N50'

  input:
  file ('*')

  
  """
  abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --out_prefix salmon --name_sample_by_basedir salmon/quant.sf
  """ 
  }

process CONTIG_EXN50_STATISTIC {

  input:

  output:
  
  """
  contig_ExN50_statistic.pl salmon.isoform.TPM.not_cross_norm {input_fasta} | tee ExN50.stats
  """ 
  }


//*Paired end rule
process ALIGN_AND_ESTIMATE_PAIRED {

  input:

  output:
  
  """
  abundance_dir = f'salmon'
  align_and_estimate_abundance.pl --thread_count {threads} --transcripts {input_fasta} --seqType fq --left {leftReads} --right {rightReads} --est_method salmon --trinity_mode --output_dir {abundance_dir}
  """ 
  }
//
    