// DATABASES

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use.
It also comes with a "auto-download" if a database is not available. Doing it the following way:
1. take userinput DB 2. add a cloud preload DB if cloud profile 3. check if the preload file exists (local or cloud)
4. if nothing is true -> download the DB and store it in the "preload" section (either cloud or local for step 3.)
*/


process get_database {
  storeDir "nextflow-autodownload-databases/busco/${params.busco}"
  stageOutMode 'move'
  containerOptions = "--user root"
  label 'BASICS'

  output:
    file("${params.busco}.tar.gz")
    
  script:
    """
    wget --quiet http://busco.ezlab.org/v2/datasets/${params.busco}.tar.gz 
    #tar -xzf ${params.busco}.tar.gz 
    #rm ${params.busco}.tar.gz
    #rm ${params.busco}/${params.busco}.tar.gz
    """
}