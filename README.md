<!--[![Build Status](https://travis-ci.org/hoelzer/ribap.svg?branch=master)](https://travis-ci.org/hoelzer/ribap)-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Nextflow](https://img.shields.io/badge/Language-Nextflow-green.svg)
[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/martinhoelzer?label=%40martinhoelzer&style=social)](https://twitter.com/martinhoelzer)

# TREAT: TranscRiptome EvaluATion

_Treat your assemblies well!_

A pipeline combining docker containers or conda environments and nextflow for the evaluation of (_de
novo_) transcriptome assembly results.

## Execution

To see all options:
````
nextflow run main.nf --help
````

Simple example execution:
````
nextflow run main.nf --assemblies test_data/rna-spades.fasta --reads test_data/eco.fastq --threads 2 --busco bacteria_odb9 -profile standard
````

## Profiles

The workflow can be deployed in different environments. 

* __standard__: uses docker containers and can be run locally
* __conda__: uses conda environments and can be run locally
* __googlegenomics__: uses the google cloud and needs to configured manually to run with your gcloud environment


## Motivation

This workflow implementation was motivated by our study comparing various _de novo_ transcriptome assembly tools:

* [Hölzer, Martin, and Manja Marz. "De novo transcriptome assembly: A comprehensive cross-species comparison of short-read RNA-Seq assemblers." GigaScience 8.5 (2019): giz039.](https://doi.org/10.1093/gigascience/giz039)
* [Online resource](https://www.rna.uni-jena.de/supplements/assembly/)