/*
* detonate call
*/
process DETONATE {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  '''
  rsem-prepare-reference --bowtie {args.refTranscripts} {args.outDir}/rsem_ref
  ''' 
}

process DETONATE {
  label 'DETONATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  """
  rsem-prepare-reference --bowtie {args.refTranscripts} {args.outDir}/rsem_ref
  """
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
  rsem-prepare-reference --bowtie {args.refTranscripts} {args.outDir}/rsem_ref
  """
}

/* Comments:
rnaQUAST is still on python2 in the rnaQUAST.py executable the shebang is still which defaults to #!/usr/bin/env python, but should now default to #!/usr/bin/env python2, because python3 is default on new systems...
*/


logger.info('SINGLE END ASSEMBLY')

    ##########################
    ## ESTIMATE THE TRUE ASSEMBLY
    # construct estimate of "true" assembly with RSEM
    logger.info('RSEM-PREPARE-REFERENCE')
    os.system(f'rsem-prepare-reference --bowtie {args.refTranscripts} {args.outDir}/rsem_ref')
    logger.info('RSEM-CALCULATE-EXPRESSION')
    os.system('rsem-calculate-expression -p {args.threads} {args.reads} {args.outDir}/rsem_ref {args.outDir}/rsem_expr')

    # command for Samtools 0.18 from the conda detonate package doesnt support -T and --threads
    logger.info('SAMTOOLS SORT')
    os.system(f'samtools sort {args.outDir}/rsem_expr.transcript.bam {args.outDir}/rsem_expr.transcript.sorted')

    # Now we are ready to estimate the “true” assembly:
    logger.info('REF-EVAL-ESTIMATE-TRUE-ASSEMBLY')
    os.system(f'ref-eval-estimate-true-assembly --reference {args.outDir}/rsem_ref --expression {args.outDir}/rsem_expr --assembly {args.outDir}/ta --alignment-policy best')

    ##########################
    ## COMPUTE KMER-COMPRESSION SCORE FOR EACH ASSEMBLY

    # This time, we will run RSEM to estimate the expression levels of each sequence in the estimated “true” assembly, as follows
    logger.info('RSEM-PREPARE-REFERENCE-2')
    os.system(f'rsem-prepare-reference --bowtie {args.outDir}/ta_0.fa {args.outDir}/ta_0_ref')
    logger.info('RSEM-CALCULATE-EXPRESSION-2')
    os.system(f'rsem-calculate-expression -p {args.threads} {args.reads} {args.outDir}/ta_0_ref {args.outDir}/ta_0_expr')

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

    for assembly in assemblies:
        assembly = abspath(assembly)  # absolute path
        tool = basename_without_ext(assembly)  # basename without ext

        logger.info(f'RSEM-EVAL-CALCULATE-SCORE: {tool}')
        os.system(f'rsem-eval-calculate-score -p {args.threads} --transcript-length-parameters {args.outDir}/transcript_length_parameters.txt {args.reads} {assembly} {args.outDir}/rsem_eval_{tool} {readLength}')

        # We now compute the KC score of each assembly as follows:
        logger.info(f'REF-EVAL-KC (KC score): {tool}')

        output = os.popen(f'ref-eval --scores kc --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --B-expr {args.outDir}/ta_0_expr.isoforms.results --kmerlen {args.readLength} --readlen {args.readLength} --num-reads {numReads}')

        kcFile = f'{args.outDir}/kc_{tool}.txt'
        with open(kcFile, "w+") as writer:
            writer.write(output.read())

        # COMPUTE ALIGNMENT-BASED SCORES FOR EACH ASSEMBLY WITH BLAT
        logger.info(f'BLAT (Alignment based scores): {tool}')
        os.system(f'blat -minIdentity=80 {args.outDir}/ta_0.fa {assembly} {args.outDir}/{tool}_to_ta_0.psl')
        os.system(f'blat -minIdentity=80 {assembly} {args.outDir}/ta_0.fa {args.outDir}/ta_0_to_{tool}.psl')

        # We can now compute the contig and nucleotide scores as follows:
        logger.info(f'REF-EVAL (Contig/Nucleotide scores): {tool}')
        output = os.popen(f'ref-eval --scores contig,nucl --weighted no --A-seqs {assembly} --B-seqs {args.outDir}/ta_0.fa --A-to-B {args.outDir}/{tool}_to_ta_0.psl --B-to-A {args.outDir}/ta_0_to_{tool}.psl --min-frac-identity 0.90')

        contigNuclFile = f'{args.outDir}/contig_nucl_{tool}.txt'
        with open(contigNuclFile, "w+") as writer:
            writer.write(output.read())