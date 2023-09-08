# Dockerized Genomic Processing Pipelines in WDL
Genomic data processing pipelines written in WDL making use of Docker containers to ensure stability. The Docker containers used are curated for usage in the Faryabi Lab processing pipelines, and can be found on [Dockerhub](https://hub.docker.com/u/faryabilab).

## Setup
#### Dependencies
Ensure that the latest `cromwell` engine is installed from their [GitHub repository](https://github.com/broadinstitute/cromwell). \
Have [Docker](https://www.docker.com/products/personal/) installed on your machine.
#### Workflows
`git clone` this repository to your local machine. **Be aware** that workflows use relative imports for tasks and sample-level workflows that are dependent on the directory structure of this repository. If you choose to move a workflow and/or task files to another location, please make sure to change the imports within the workflow source files directly, if needed.
#### Sample Sheet
Each run of a workflow must be accompanied by a `samplesheet.tsv` - a single-column tab-separated file that specifies the prefixes ('sample names') of each `fastq` file to be processed. \
**A note on fastq naming** - Ensure that your paired-end `fastq` files follow the naming convention `XYZ_R1.fastq.gz`, where `XYZ` represents this file's sample name. Single-end `fastq`s must be named as `XYZ.fastq.gz`.
To make a sample sheet, simply run `make_samplesheet.sh` found in `utils/`: 
* Paired-end: `./make_samplesheet.sh -p -d /path/to/fastq_dir`
* Single-end: `./make_samplesheet.sh -d /path/to/fastq_dir` \

This will create a `samplesheet.tsv` in your current working directory.
## Cromwell Configuration
Within the `cromwell_configs` directory is a Cromwell config file which instructs the engine to use Docker as backend. Right now, **the only option you should change** is `concurrent-job-limit`.
* `concurrent-job-limit` - Default = 10, most amount of jobs that can be run at once.
#### Call-Caching
A desirable feature for many workflows is the ability to identify steps, or "jobs", that have already been run and essentially "skip" that step, saving both time and computing resources. In this reporisoty, Cromwell is configured to use a local MySQL database (which runs from a Docker container) called `PipelineDatabase`. If you want to use a different SQL, or any, database, the `database` block must be edited accordingly, else, the config will attempt to log in to the MySQL server specified. 
## Usage
First, navigate to the workflow you want to use, located the in `./workflows` directory. Each subdirectory contains the batch and sample-level workflows, as well as the config file. The `*_inputs.config` file must be configured for the workflow to run properly. Many of these options will be parameters for the command-line tools used in the workflow, but the important ones are: 
* `project_out_dir` - Main output directory, where by-sample directories will be created for by-cample output files.
* `fastq_dir` - Path to one directiry with **all** Fastq files to be processed.
* `sampleList` - Samplesheet with one column of sample names and no header. Can be created with utility script `utils/make_samplesheet.sh`. See Setup section for details.
* `paired` - Boolean, `true` for paired-end experiemnts, `false` otherwise.
* `ChromNoScaffold` - BED-style file with enire chromosome intervals to keep in the resulting BAM.
* `ChromosomeSizes` - BEDtools-style chromosome sizes file.
* `Blacklist` - A BED-style file containing regions to be removed from the resulting BAM, typically known problematic regions.
Once the configuration is complete, one can run a workflow as such: 
```
java -Dconfig.file=/path/to/cromwell_config -jar /path/to/cromwell.jar run batch_workflow.wdl -i workflow_input.json
```
## Scalability 
Currently, the pipeline's scalable features are implemented as `cpu` and `mem` input options in `<workflow>_inputs.json`. Each memory/compute-heavy tasks is given its own version of the variable, with `cpu` representing cores and `mem` representing memory in gigabytes. \
The aforementioned `concurrent-job-limit` variable within the Cromwell config file sets the upper limit on how many jobs, or samples, that can be processed at once.






