/*
* Benchmark universal orthologous groups w/ BUSCO
*/ 
process BUSCO {
  label 'BUSCO'
  publishDir "${params.output}/${params.dir}/", mode: 'copy', pattern: "${name}_busco_summary.txt"
  publishDir "${params.output}/${params.dir}/", mode: 'copy', pattern: "${name}_busco_figure.pdf"

  input:
  set val(name), file(assembly)
  file(db)

  output:
  set val(name), file("${name}_busco_summary.txt"), file ("${name}_busco_figure.pdf")

  script:
  """
  export AUGUSTUS_CONFIG_PATH="/opt/conda/config/"
  # run BUSCO
  run_BUSCO.py -i ${assembly} -o results -l ${db} -m tran -c ${params.threads} -t ./ -z
  cp run_results/short_summary_results.txt ${name}_busco_summary.txt
  # generate Plot and rehack Rscript
  generate_plot.py -wd run_results/
  sed -i 's/busco_figure.png/busco_figure.pdf/g' run_results/busco_figure.R
  Rscript run_results/busco_figure.R
  cp run_results/busco_figure.pdf  ${name}_busco_figure.pdf
  """
}

/* Comments:
Buscos plotting feature does not work in Docker by default.
The "hack" modifies via sed the Rscript that gets generated after generate_plot.py
After changing to pdf (circumventing resolutions) we rerun the script via Rscript. Tadaa... we have a plot
*/
