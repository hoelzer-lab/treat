process buscoGetDB {
  storeDir 'nextflow-autodownload-databases/busco'
  stageOutMode 'move'
  label 'BUSCO'
  containerOptions = "--user root"

  
  input:
    val(db)

  output:
    file(db)
    
  script:
    """
    wget --quiet http://busco.ezlab.org/v2/datasets/${db}.tar.gz 
    tar -xzf ${db}.tar.gz 
    rm ${db}.tar.gz
    """
}