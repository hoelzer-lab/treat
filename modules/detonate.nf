/*
* DETONATE FOR SINGLE END READS
*/

// args.outDir = params.dir

workflow DETONATE {

  main:
    // performed once   
    ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS(transcripts_ch)
    RSEM_PREPARE_REFERENCE(transcripts_ch)
    RSEM_CALCULATE_EXPRESSION(RSEM_PREPARE_REFERENCE.out, reads_detonate)
    SAMTOOLS_SORT(RSEM_CALCULATE_EXPRESSION.out.filter{ "~/*.bam/" })
    REF_EVAL_ESTIMATE_TRUE_ASSEMBLY(RSEM_PREPARE_REFERENCE.out, RSEM_CALCULATE_EXPRESSION.out, SAMTOOLS_SORT.out)
    RSEM_PREPARE_REFERENCE_2(REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    RSEM_CALCULATE_EXPRESSION_2(RSEM_PREPARE_REFERENCE_2.out, reads_detonate)
    DETONATE7(reads_detonate, SAMTOOLS_SORT.out)
    // performed for each assembly
    RSEM_EVAL_CALCULATE_SCORE(ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS.out, reads_detonate, assemblies_ch)
    REF_EVAL_KC(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out, RSEM_CALCULATE_EXPRESSION_2.out)
    BLAT_A_TO_B(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    BLAT_B_TO_A(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    REF_EVAL_CONTIG(assemblies_ch, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out, BLAT_A_TO_B.out, BLAT_B_TO_A.out)
  emit:
    RSEM_CALCULATE_EXPRESSION.out
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

process DETONATE7 {
  label 'DETONATE'

  input:
  file(reads)
  file(sorted_bam)

  output:
  stdout()
  
  shell:
  """
    # Estimate the num-reads parameter, this is the number of reads times the read length equals the total number of nucleotides in the read set.
    echo awk '{s++}END{print s/4}' ${reads}

    # Get read length and insert size info via the mapped bam file from bowtie
    #samtools stats rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' | cut -f 3
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

  
  output:
  file("rsem_eval*")
  
  shell:
  """
  rsem-eval-calculate-score -p !{params.threads} --transcript-length-parameters transcript_length_parameters.txt ${reads} ${assembly} rsem_eval_${name} 100
  #rsem-eval-calculate-score -p !{params.threads} --transcript-length-parameters transcript_length_parameters.txt ${reads} {assembly} {args.outDir}/rsem_eval_{tool} {readLength}
  """
}

process REF_EVAL_KC {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "kc_${name}.txt"

  input:
  tuple val(name), file(assembly)
  file('ta_0.fa')
  file(test)
  
  output:
  file("kc_${name}.txt")
  
  shell:
  """
  ref-eval --scores kc --A-seqs ${assembly} --B-seqs ta_0.fa --B-expr ta_0_expr.isoforms.results --kmerlen 75 --readlen 100 --num-reads 100 | tee kc_${name}.txt
  #ref-eval --scores kc --A-seqs ${assembly} --B-seqs {args.outDir}/ta_0.fa --B-expr {args.outDir}/ta_0_expr.isoforms.results --kmerlen {args.readLength} --readlen {args.readLength} --num-reads {numReads} | tee kc_${name}.txt
  """
}

process BLAT_A_TO_B {
  label 'DETONATE'

  input:
  tuple val(name), file(assembly)
  file('ta_0.fa')

  output:
  file("${name}_to_ta_0.psl")
  
  shell:
  """
  blat -minIdentity=80 ta_0.fa ${assembly} ${name}_to_ta_0.psl
  """
}

process BLAT_B_TO_A {
  label 'DETONATE'

  input:
  tuple val(name), file(assembly)
  file('ta_0.fa')

  output:
  file("ta_0_to_${name}.psl")
  
  shell:
  """
  blat -minIdentity=80 ${assembly} ta_0.fa ta_0_to_${name}.psl
  """
}

process REF_EVAL_CONTIG {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "contig_nucl_${name}.txt"

  input:
  tuple val(name), file(assembly)
  file('ta_0.fa')
  file(test)
  file(test2)
  

  output:
  file("contig_nucl_${name}.txt")
  
  shell:
  """
  ref-eval --scores contig,nucl --weighted no --A-seqs ${assembly} --B-seqs ta_0.fa --A-to-B ${name}_to_ta_0.psl --B-to-A /ta_0_to_${name}.psl --min-frac-identity 0.90 | tee contig_nucl_${name}.txt
  """
}

// /* Comments:
// */