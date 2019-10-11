process buscoGetDB {
  storeDir 'nextflow-autodownload-databases/busco'
  label 'BUSCO'
  
  input:
    val(db)

  output:
    file(db)
    
  script:
    """
    wget http://busco.ezlab.org/v2/datasets/${db}.tar.gz 
    tar -xvzf ${db}.tar.gz
    rm ${db}.tar.gz
    """
}