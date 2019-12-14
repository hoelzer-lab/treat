/*
* Benchmark universal orthologous groups w/ BUSCO
*/ 

// get BUSCO db
include './buscoGetDB' params(busco: params.busco)

workflow BUSCO{
    main:
        get_database()
        run_busco(params.assemblies, get_database.out)
        plot_busco(run_busco.out.collect())
    emit:
        run_busco.out
}

process run_busco {
  label 'BUSCO'
  publishDir "${params.output}/${params.dir}/", mode: 'copy', pattern: "short_summary_${name}.txt"
  // publishDir "${params.output}/${params.dir}/", mode: 'copy', pattern: "${name}_busco_figure.pdf"

  input:
  tuple val(name), file(assembly)
  file(db)

  output:
  file("short_summary_${name}.txt")

  shell:
  """
  tar -xzf ${db}
  export AUGUSTUS_CONFIG_PATH="/opt/conda/config/"
  # run BUSCO
  run_BUSCO.py -i ${assembly} -o ${name} -l ${db.simpleName} -m tran -c ${params.threads} -t ./ -z
  cp run_${name}/short_summary_${name}.txt short_summary_${name}.txt
  """
}

process plot_busco {
  label 'BUSCO'
  publishDir "${params.output}/${params.dir}/", mode: 'copy', pattern: "busco_figure.pdf"

  input:
  file('*')

  output:
  file ("busco_figure.pdf")

  script:
  """
  export AUGUSTUS_CONFIG_PATH="/opt/conda/config/"
  mkdir run_results
  cp *.txt run_results/
  # generate Plot and rehack Rscript
  generate_plot.py -wd run_results/
  sed -i 's/busco_figure.png/busco_figure.pdf/g' run_results/busco_figure.R
  Rscript run_results/busco_figure.R
  cp run_results/busco_figure.pdf busco_figure.pdf
  """
}

/* Comments:
Buscos plotting feature does not work in Docker by default.
The "hack" modifies via sed the Rscript that gets generated after generate_plot.py
After changing to pdf (circumventing resolutions) we rerun the script via Rscript. Tadaa... we have a plot
*/
