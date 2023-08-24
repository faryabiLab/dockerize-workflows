# Dockerized Genomic Processing Pipelines in WDL
Genomic data processing pipelines written in WDL making use of Docker containers to ensure stability.

## Setup
`git clone` the repository to your local machine. Ensure that the latest `cromwell` engine is installed from their [GitHub repository](https://github.com/broadinstitute/cromwell).
## Usage
First, navigate to the workflow you want to use, located the in `./workflows` directory. Each subdirectory contains the batch and sample-level workflows, as well as the config file. The `*_inputs.config` file must be configured for the workflow to run properly. Many of these options will be parameters for the command-line tools used in the workflow, but the important ones are: 
* `project_out_dir` - Main output directory, where by-sample directories will be created for by-cample output files.
* `fastq_dir` - Path to one directiry with **all** Fastq files to be processed.
* `sampleList` - Samplesheet with one column of sample names. Can be created with utility script `utils/make_samplesheet.sh`.

Once the `inputs.json` file has been properly configured, the workflow can be run as such:
```
java -Dconfig.file=cromwell_configs/cromwell.docker.config -jar cromwell-xx.jar run workdlow.wdl -i inputs.json
```
