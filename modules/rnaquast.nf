/*
* rnaquast call
*/
process RNAQUAST {
  label 'RNAQUAST'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  '''
  rnaQUAST.py --help
  ''' 
}

/* Comments:
rnaQUAST is still on python2 in the rnaQUAST.py executable the shebang is still which defaults to #!/usr/bin/env python, but should now default to #!/usr/bin/env python2, because python3 is default on new systems...
*/
