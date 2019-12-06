/*
* DETONATE FOR SINGLE END READS
*/

// args.outDir = params.dir

workflow DETONATE {

  main:
    // performed once   
    ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS(params.transcripts)
    RSEM_PREPARE_REFERENCE(params.transcripts)
    RSEM_CALCULATE_EXPRESSION(RSEM_PREPARE_REFERENCE.out, params.reads)
    SAMTOOLS_SORT(RSEM_CALCULATE_EXPRESSION.out.filter{ "~/*.bam/" })
    REF_EVAL_ESTIMATE_TRUE_ASSEMBLY(RSEM_PREPARE_REFERENCE.out, RSEM_CALCULATE_EXPRESSION.out, SAMTOOLS_SORT.out)
    RSEM_PREPARE_REFERENCE_2(REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    RSEM_CALCULATE_EXPRESSION_2(RSEM_PREPARE_REFERENCE_2.out, params.reads)
    NUMBER_OF_READS(params.reads)
    READ_LENGTH(SAMTOOLS_SORT.out)
    // performed for each assembly
    RSEM_EVAL_CALCULATE_SCORE(ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS.out, params.reads, params.assemblies, READ_LENGTH.out)
    REF_EVAL_KC(params.assemblies, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out, RSEM_CALCULATE_EXPRESSION_2.out, NUMBER_OF_READS.out, READ_LENGTH.out)
    //BLAT_A_TO_B(params.assemblies, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    //BLAT_B_TO_A(params.assemblies, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    BLAT(params.assemblies, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    REF_EVAL_CONTIG(params.assemblies, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out, BLAT.out)
  emit:
    kc = REF_EVAL_KC.out
    contig = REF_EVAL_CONTIG.out
    rsem = RSEM_EVAL_CALCULATE_SCORE.out
}



process ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS {
  label 'DETONATE'

  input:
  file(transcripts)

  output:
  file('transcript_length_parameters.txt')

  shell:
  """
  rsem-eval-estimate-transcript-length-distribution ${transcripts} transcript_length_parameters.txt
  """
}

process RSEM_PREPARE_REFERENCE {
  label 'DETONATE'

  input:
  file(transcripts)

  output:
  file("rsem_ref*")
  
  shell:
  """
  rsem-prepare-reference --bowtie ${transcripts} rsem_ref
  """ 
}

process RSEM_CALCULATE_EXPRESSION {
  label 'DETONATE'

  input:
  file(assemblies)
  file(reads)
  
  output:
    file("rsem_expr*")
  
  shell:
  """
  rsem-calculate-expression -p !{params.threads} ${reads} rsem_ref rsem_expr
  """
}

process SAMTOOLS_SORT {
  label 'DETONATE'

  input:
  file(sorted_bam)
  
  output:
  file("rsem_expr.transcript.sorted.bam")
  
  shell:
  """
  samtools sort --threads !{params.threads} -o rsem_expr.transcript.sorted.bam rsem_expr.transcript.bam
  """
}

process REF_EVAL_ESTIMATE_TRUE_ASSEMBLY {
  label 'DETONATE'

  input:
  file(rsem_ref)
  file(rsem_expr)
  file(sorted_bam)

  output:
  file('ta_0.fa')
  
  shell:
  """
  ref-eval-estimate-true-assembly --reference rsem_ref --expression rsem_expr --assembly ta --alignment-policy best
  """
}

process RSEM_PREPARE_REFERENCE_2 {
  label 'DETONATE'

  input:
  file('ta_0.fa')  

  output:
  file('ta_0_ref*')
  
  shell:
  """
  rsem-prepare-reference --bowtie ta_0.fa ta_0_ref
  """
}

process RSEM_CALCULATE_EXPRESSION_2 {
  label 'DETONATE'

  input:
  file('*')
  file(reads)
  
  output:
  file('ta_0_expr*')
  
  shell:
  """
  rsem-calculate-expression -p !{params.threads} ${reads} ta_0_ref ta_0_expr
  """
}

process NUMBER_OF_READS {
  label 'DETONATE'

  input:
  file(reads)

  output:
  stdout()
  
  shell:
  """
    # Estimate the num-reads parameter, this is the number of reads times the read length equals the total number of nucleotides in the read set.
    awk '{s++}END{printf s/4}' ${reads}
  """
}

process READ_LENGTH {
  label 'DETONATE'

  input:
  file(sorted_bam)

  output:
  stdout()
  
  shell:
  """
    # Get read length and insert size info via the mapped bam file from bowtie
    samtools stats rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' | awk -F" " '{printf(\$4)}'
  """
}

// THESE METHODS ARE FOR EACH ASSEMBLY

process RSEM_EVAL_CALCULATE_SCORE {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "rsem_eval_${name}.score"

  input:
  file('transcript_length_parameters.txt')
  file(reads)
  tuple val(name), file(assembly)
  val (read_length)

  output:
    file("rsem_eval*")
  
  shell:
  """
  rsem-eval-calculate-score -p !{params.threads} --transcript-length-parameters transcript_length_parameters.txt ${reads} ${assembly} rsem_eval_${name} ${read_length}
  """
}

process REF_EVAL_KC {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "kc_${name}.txt"

  input:
  tuple val(name), file(assembly)
  file('*')
  file(test)
  val (number_of_reads)
  val (read_length)
  
  output:
  file("kc_${name}.txt")
  
  shell:
  """
  ref-eval --scores kc --A-seqs ${assembly} --B-seqs ta_0.fa --B-expr ta_0_expr.isoforms.results --kmerlen ${read_length} --readlen ${read_length} --num-reads ${number_of_reads} | tee kc_${name}.txt
  """
}

process BLAT {
  label 'DETONATE'

  input:
  tuple val(name), file(assembly)
  file('*')

  output:
  tuple file("ta_0_to_${name}.psl"), file("${name}_to_ta_0.psl")
  
  shell:
  """
  blat -minIdentity=80 ${assembly} ta_0.fa ta_0_to_${name}.psl
  blat -minIdentity=80 ta_0.fa ${assembly} ${name}_to_ta_0.psl
  """
}

process REF_EVAL_CONTIG {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "contig_nucl_${name}.txt"

  input:
  tuple val(name), file(assembly)
  file('*')
  file("*")

  output:
  file("contig_nucl_${name}.txt")
  
  shell:
  """
  ref-eval --scores contig,nucl --weighted no --A-seqs ${assembly} --B-seqs ta_0.fa --A-to-B ${name}_to_ta_0.psl --B-to-A ta_0_to_${name}.psl --min-frac-identity 0.90 > contig_nucl_${name}.txt
  """
}

// /* Comments:
// */