/*
* Transrate call
*/
process TRANSRATE {
  label 'TRANSRATE'
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/assemblies.csv"

  input:
  set val(name), file(assembly)

  output:
  set val(name), file("${name}/assemblies.csv")
  
  shell:
  '''
  transrate --output !{name} --assembly !{assembly}
  ''' 
}

/* Comments:
*/
