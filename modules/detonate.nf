/*
* DETONATE FOR SINGLE END READS
*/

// args.outDir = params.dir
process DETONATE {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  '''
  rsem-prepare-reference --bowtie {args.refTranscripts} !{params.dir}/rsem_ref
  ''' 
}

process DETONATE2 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  rsem-calculate-expression -p !{params.threads} {args.reads} {args.outDir}/rsem_ref {args.outDir}/rsem_expr
  """
}

process DETONATE3 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
samtools sort {args.outDir}/rsem_expr.transcript.bam {args.outDir}/rsem_expr.transcript.sorted
  """
}

process DETONATE4 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  ref-eval-estimate-true-assembly --reference {args.outDir}/rsem_ref --expression {args.outDir}/rsem_expr --assembly {args.outDir}/ta --alignment-policy best
  """
}

process DETONATE5 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  rsem-prepare-reference --bowtie {args.outDir}/ta_0.fa {args.outDir}/ta_0_ref
  """
}

process DETONATE6 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  rsem-calculate-expression -p !{params.threads} {args.reads} {args.outDir}/ta_0_ref {args.outDir}/ta_0_expr
  """
}


process DETONATE7 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
    # Estimate the num-reads parameter, this is the number of reads times the read length equals the total number of nucleotides in the read set.
    output = os.popen("awk '{s++}END{print s/4}' " + args.leftReads)
    numReads = str(output.read())

    # Get read length and insert size info via the mapped bam file from bowtie
    p = subprocess.run(f"samtools stats {args.outDir}/rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' | cut -f 3", stdout=subprocess.PIPE, shell=True)
    readLength = str(p.stdout.decode('utf-8'))
    logger.info(f"Read length: {readLength}")
  """
}


process DETONATE8 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
    # Estimate the num-reads parameter, this is the number of reads times the read length equals the total number of nucleotides in the read set.
    logger.info('Calculating number of reads.')
    output = os.popen("awk '{s++}END{print s/4}' " + args.leftReads)
    numReads = str(output.read())
    logger.info(f"Number of reads: {numReads}")

    # Get read length and insert size info via the mapped bam file from bowtie
    logger.info('Calculating read length and insert size.')
    p = subprocess.run(f"samtools stats {args.outDir}/rsem_expr.transcript.sorted.bam | grep -e '^SN\tmaximum length' | cut -f 3", stdout=subprocess.PIPE, shell=True)
    readLength = str(p.stdout.decode('utf-8'))
    logger.info(f"Read length: {readLength}")
  """
}

// FOR EACH ASSEMBLY PROVIDED:
//assembly = abspath(assembly)  # absolute path
        //tool = basename_without_ext(assembly)  # basename without ext

process DETONATE9 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  rsem-eval-calculate-score -p !{params.threads} --transcript-length-parameters {args.outDir}/transcript_length_parameters.txt {args.reads} {assembly} {args.outDir}/rsem_eval_{tool} {readLength}
  """
}

process DETONATE10 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  output = os.popen(f'ref-eval --scores kc --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --B-expr {args.outDir}/ta_0_expr.isoforms.results --kmerlen {args.readLength} --readlen {args.readLength} --num-reads {numReads}')

  kcFile = f'{args.outDir}/kc_{tool}.txt'
  with open(kcFile, "w+") as writer:
      writer.write(output.read())
  """
}


process DETONATE11 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  blat -minIdentity=80 {args.outDir}/ta_0.fa {assembly} {args.outDir}/{tool}_to_ta_0.psl
  """
}


process DETONATE12 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  blat -minIdentity=80 {assembly} {args.outDir}/ta_0.fa {args.outDir}/ta_0_to_{tool}.psl
  """
}

process DETONATE13 {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)
  

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  output = os.popen(f'ref-eval --scores contig,nucl --weighted no --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --A-to-B {args.outDir}/{tool}_to_ta_0.psl --B-to-A {args.outDir}/ta_0_to_{tool}.psl --min-frac-identity 0.90')

  contigNuclFile = f'{args.outDir}/contig_nucl_{tool}.txt'
  with open(contigNuclFile, "w+") as writer:
    writer.write(output.read())
  """
}

/* Comments:
rnaQUAST is still on python2 in the rnaQUAST.py executable the shebang is still which defaults to #!/usr/bin/env python, but should now default to #!/usr/bin/env python2, because python3 is default on new systems...
*/