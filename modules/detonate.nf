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
    BLAT(params.assemblies, REF_EVAL_ESTIMATE_TRUE_ASSEMBLY.out)
    REF_EVAL_CONTIG(BLAT.out)
  emit:
    kc = REF_EVAL_KC.out
    contig = REF_EVAL_CONTIG.out
    rsem = RSEM_EVAL_CALCULATE_SCORE.out
}



process ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS {
  label "DETONATE"

  input:
    file(transcripts)

  output:
    file("transcript_length_parameters.txt")

  shell:
    """
    rsem-eval-estimate-transcript-length-distribution ${transcripts} transcript_length_parameters.txt
    """
}

process RSEM_PREPARE_REFERENCE {
  label "RSEM"

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
  label "RSEM"

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
  label "DETONATE"

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
  label "DETONATE"

  input:
    file(rsem_ref)
    file(rsem_expr)
    file(sorted_bam)

  output:
    file("ta_0.fa")
  
  shell:
    """
    ref-eval-estimate-true-assembly --reference rsem_ref --expression rsem_expr --assembly ta --alignment-policy best
    """
}

process RSEM_PREPARE_REFERENCE_2 {
  label "RSEM"

  input:
    file("ta_0.fa")

  output:
    file("ta_0_ref*")
  
  shell:
    """
    rsem-prepare-reference --bowtie ta_0.fa ta_0_ref
    """
}

process RSEM_CALCULATE_EXPRESSION_2 {
  label "RSEM"

  input:
    file("*")
    file(reads)
  
  output:
    file("ta_0_expr*")
  
  shell:
    """
    rsem-calculate-expression -p !{params.threads} ${reads} ta_0_ref ta_0_expr
    """
}

process NUMBER_OF_READS {
  label "DETONATE"

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
  label "DETONATE"

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
  label "DETONATE"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "rsem_eval_${name}.score"

  input:
    file("transcript_length_parameters.txt")
    file(reads)
    tuple val(name), file(assembly)
    val (read_length)

  output:
    file("rsem_eval*")
  
  shell:
    """
    rsem-eval-calculate-score ${reads} ${assembly} rsem_eval_${name} ${read_length} --transcript-length-parameters transcript_length_parameters.txt -p !{params.threads}
    """
}

process REF_EVAL_KC {
  label "DETONATE"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "kc_${name}.txt"

  input:
    tuple val(name), file(assembly)
    file("*")
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
  label "DETONATE"

  input:
    tuple val(name), file(assembly)
    file("ta_0.fa")

  output:
    tuple file("ta_0_to_${name}.psl"), file("${name}_to_ta_0.psl")
    tuple val(name), file(assembly)
    file("ta_0.fa")
  
  shell:
    """
      blat -minIdentity=80 ${assembly} ta_0.fa ta_0_to_${name}.psl
      blat -minIdentity=80 ta_0.fa ${assembly} ${name}_to_ta_0.psl
    """
}

process REF_EVAL_CONTIG {
  label "DETONATE"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "contig_nucl_${name}.txt"

  input:
    tuple file(ta_to_assembly), file(assembly_to_ta)
    tuple val(name), file(assembly)
    file("ta_0.fa")


  output:
    file("contig_nucl_${name}.txt")
  
  shell:
    """
      ref-eval --scores contig,nucl --weighted no --A-seqs ${assembly} --B-seqs ta_0.fa --A-to-B ${assembly_to_ta} --B-to-A ${ta_to_assembly} --min-frac-identity 0.90 > contig_nucl_${name}.txt
    """

}

/*
##################################################
########## DETONATE PAIRED END PIPELINE ##########
##################################################
The process ESTIMATE_TRANSCRIPTS_LENGTH_PARAMETERS is needed at the beginning of both paired- and single-end pipelines.

*/
process PE_RSEM_PREPARE_REFERENCE {
  label "RSEM"

  input:
    file(transcripts)

  output:
    file("rsem_ref*")
  
  shell:
    """
    rsem-prepare-reference --bowtie ${transcripts} rsem_ref
    """ 
}

process PE_RSEM_CALCULATE_EXPRESSION {
  label "RSEM"

  input:
    file(assemblies)
    tuple file(leftReads), file(rightReads)
  
  output:
    file("rsem_expr*")
  
  shell:
    """
    rsem-calculate-expression -p !{params.threads} --paired-end ${leftReads} ${rightReads} rsem_ref rsem_expr")
    """
}

process PE_SAMTOOLS_SORT {
  label "DETONATE"

  input:
    file(sorted_bam)
  
  output:
    file("rsem_expr.transcript.sorted.bam")
  
  shell:
    """
    samtools sort --threads !{params.threads} -o rsem_expr.transcript.sorted.bam rsem_expr.transcript.bam
    """
}

process PE_REF_EVAL_ESTIMATE_TRUE_ASSEMBLY {
  label "DETONATE"

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

process PE_RSEM_PREPARE_REFERENCE_2 {
  label "RSEM"

  input:
    file("ta_0.fa")

  output:
    file("ta_0_ref*")
  
  shell:
    """
    rsem-prepare-reference --bowtie ta_0.fa ta_0_ref
    """
}

process PE_RSEM_CALCULATE_EXPRESSION_2 {
  label "RSEM"

  input:
    file("*")
    tuple file(leftReads), file(rightReads)
  
  output:
    file("ta_0_expr*")
  
  shell:
    """
    rsem-calculate-expression -p !{params.threads} --paired-end ${leftReads} ${rightReads} ta_0_ref ta_0_expr
    """
}

process PE_NUMBER_OF_READS {
  label "DETONATE"

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

process PE_READ_LENGTH_AND_INSERT_SIZE {
  // TODO: Create a little python script that calculates the fragment length from the insert size and the read length
  // fragment size = 2xreadLength + insertSize
  // p = subprocess.run(f"samtools stats {args.outDir}/rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' -e '^SN\tinsert size average' | cut -f 3", stdout=subprocess.PIPE, shell=True)
  // output = p.stdout.decode('utf-8').split()
  // readLength = int(round(float(output[0])))
  // insertSize = float(output[1])
  // # for paired end reads: fragment size = 2xreadLength + insertSize
  // fragmentLength = (2 * readLength) + insertSize
  label "DETONATE"

  input:
    file(sorted_bam)

  output:
    stdout()
  
  shell:
  """
    # Get read length and insert size info via the mapped bam file from bowtie
    samtools stats rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' -e '^SN\tinsert size average' | cut -f 3
  """
}

// Processes for each paired-end assembly

process PE_RSEM_EVAL_CALCULATE_SCORE {
  label "DETONATE"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "rsem_eval_${name}.score"

  input:
    file('transcript_length_parameters.txt')
    tuple file(leftReads), file(rightReads)
    tuple val(name), file(assembly)
    val (read_length)

  output:
    file("rsem_eval*")
  
  shell:
    """
    rsem-eval-calculate-score --paired-end ${leftReads} ${rightReads} ${assembly} rsem_eval_${name} ${fragmentLength} --transcript-length-parameters transcript_length_parameters.txt -p !{params.threads}
    """
}
        
process PE_REF_EVAL_KC {
  label "DETONATE"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "kc_${name}.txt"

  input:
    tuple val(name), file(assembly)
    file("*")
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

process PE_BLAT {
  label "DETONATE"

  input:
    tuple val(name), file(assembly)
    file("ta_0.fa")

  output:
    tuple file("ta_0_to_${name}.psl"), file("${name}_to_ta_0.psl")
    tuple val(name), file(assembly)
    file("ta_0.fa")
  
  shell:
    """
    blat -minIdentity=80 ${assembly} ta_0.fa ta_0_to_${name}.psl
    blat -minIdentity=80 ta_0.fa ${assembly} ${name}_to_ta_0.psl
    """
}

process PE_REF_EVAL_CONTIG {
  label "DETONATE"
  publishDir "${params.output}/${params.dir}/", mode: "copy", pattern: "contig_nucl_${name}.txt"

  input:
    tuple file(ta_to_assembly), file(assembly_to_ta)
    tuple val(name), file(assembly)
    file("ta_0.fa")


  output:
    file("contig_nucl_${name}.txt")
  
  shell:
    """
    ref-eval --scores contig,nucl --weighted no --A-seqs ${assembly} --B-seqs ta_0.fa --A-to-B ${assembly_to_ta} --B-to-A ${ta_to_assembly} --min-frac-identity 0.90 > contig_nucl_${name}.txt
    """

}

 /* Comments:
 A lot of the paired-end pipeline may be redundant and can be combined with the single end pipeline. Basically an extra rule is only necessary if a mapping step is performed.
 */