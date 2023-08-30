# Dockerized Genomic Processing Pipelines in WDL
Genomic data processing pipelines written in WDL making use of Docker containers to ensure stability. The Docker containers used are curated for usage in the Faryabi Lab processing pipelines, and can be found on [Dockerhub](https://hub.docker.com/u/faryabilab).

## Setup
#### Dependencies
Ensure that the latest `cromwell` engine is installed from their [GitHub repository](https://github.com/broadinstitute/cromwell). \
Have [Docker](https://www.docker.com/products/personal/) installed on your machine. 
#### Workflows
`git clone` this repository to your local machine. **Be aware** that workflows use relative imports for tasks and sample-level workflows that are dependent on the directory structure of this repository. If you choose to move there workflow and/or task files to another location, please make sure to change the imports within the workflow source files directly, if needed.
## Usage
First, navigate to the workflow you want to use, located the in `./workflows` directory. Each subdirectory contains the batch and sample-level workflows, as well as the config file. The `*_inputs.config` file must be configured for the workflow to run properly. Many of these options will be parameters for the command-line tools used in the workflow, but the important ones are: 
* `project_out_dir` - Main output directory, where by-sample directories will be created for by-cample output files.
* `fastq_dir` - Path to one directiry with **all** Fastq files to be processed.
* `sampleList` - Samplesheet with one column of sample names. Can be created with utility script `utils/make_samplesheet.sh`.
* `paired` - Boolean, `true` for paired-end experiemnts, `false` otherwise.
* `ChromNoScaffold` - BED-style file with enire chromosome intervals to keep in the resulting BAM.
* `ChromosomeSizes` - BEDtools-style chromosome sizes fille.
* `Blacklist` - A BED-style file containing regions to be removed from the resulting BAM, typically known problematic regions.

Once the `inputs.json` file has been properly configured, the workflow can be run as such:
```
java -Dconfig.file=cromwell_configs/cromwell.docker.config -jar cromwell-xx.jar run workdlow.wdl -i inputs.json
```
