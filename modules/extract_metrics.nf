
workflow HEATMAP{
    get:
      HISAT
      BUSCO
      RNAQUAST
      DETONATE_1
      DETONATE_2
      DETONATE_3
      EX90N50
    main:
        COLLECT_INFO(HISAT, BUSCO, RNAQUAST, DETONATE_1, DETONATE_2, DETONATE_3, EX90N50)
}

process COLLECT_INFO {
  label 'HEATMAP'
  publishDir "${params.output}/", mode:'copy', pattern: "*.svg"
  publishDir "${params.output}/", mode:'copy', pattern: "*.csv"

  input:
    file(mapping)
    file(busco_metric)
    file(rnaquast)
    file(kc)
    file(contig)
    file(rsem)
    file(ex90n50)

  output:
  file("all_metrics.csv")
  file("selected_metrics.csv")
  file("selected_normalized_metrics.csv")
  file("heatmap.svg")
  
  """
  create_heatmap.py
  """ 
  }

/*
process CREATE {
  label 'HEATMAP'
  publishDir "${params.output}/", mode:'copy', pattern: "heatmap.svg"

  input:
  tuple val(name), file(mapping)
  tuple val(name), file(busco_metric), file(busco_figure)
  tuple val(name), file(rnaquast)
  file(kc)
  file(contig)
  file(rsem)

  output:
  file("heatmap.svg")
  
  """
  create_heatmap.py -m ${mapping} -b ${busco_metric} -r ${rnaquast} -name ${name}
  """ 
  }
*/