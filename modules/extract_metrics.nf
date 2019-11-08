process EXTRACT {
  publishDir "${params.output}/${params.dir}/", mode:'copy', pattern: "${name}/short_report.tsv"

  input:
  tuple val(name), file(assembly)
  tuple val(read_id), file(reads)
  file(reference)
  file(annotation)

  output:
  tuple file("heatmap.svg")
  
  shell:
  """
  """ 
}