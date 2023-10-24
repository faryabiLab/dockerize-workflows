# Dockerized Genomic Processing Pipelines in WDL
Genomic data processing pipelines written in WDL making use of Docker containers to ensure stability. The Docker containers used are curated for usage in the Faryabi Lab processing pipelines, and can be found on [Dockerhub](https://hub.docker.com/u/faryabilab).

## Setup
#### Dependencies
Ensure that the latest `cromwell` engine (.jar file) and womtools (.jar) is installed from their [GitHub repository](https://github.com/broadinstitute/cromwell). \
Have [Docker](https://www.docker.com/products/personal/) installed on your machine.
#### Sample Sheet
Each run of a workflow must be accompanied by a `samplesheet.tsv` - a single-column tab-separated file that specifies the prefixes ('sample names') of each `fastq` file to be processed. \
**A note on fastq naming** - Ensure that your paired-end `fastq` files follow the naming convention `XYZ_R1.fastq.gz` / `XYZ_R2.fastq.gz`, where `XYZ` represents this file's sample name. Single-end `fastq`s must be named as `XYZ.fastq.gz`.
To make a sample sheet, simply run `make_samplesheet.sh` found in `utils/`: 
* Paired-end: `./make_samplesheet.sh -p -d /path/to/fastq_dir`
* Single-end: `./make_samplesheet.sh -d /path/to/fastq_dir` 

This will create a `samplesheet.tsv` in your current working directory.
## Cromwell Configuration
Within the `cromwell_configs` directory is a Cromwell config file which instructs the engine to use Docker as backend. The only option that can be tweaked in this file is the `concurrent_job_limit` variable: this controls the number of jobs that can be running at once (default = 10).
#### Call-Caching
A desirable feature for many workflows is the ability to identify steps, or "jobs", that have already been run and essentially "skip" that step, saving both time and computing resources. In this reporisoty, Cromwell is configured to use a local MySQL database (which runs from a Docker container) called `PipelineDatabase`, to store intermediate files to determine the steps that have been run. If you want to use a different SQL, or any, database, the `database` block must be edited accordingly in the Cromwell configuration file, else, the config will attempt to log in to the MySQL server specified. By default, the `database` block is configured to work with a MySQL Docker container that can be initialized via the following command:
```
docker run \
-p 52000:3306 \
--name PipelineDatabase \
-e MYSQL_ROOT_PASSWORD=@noah1234 \
-e MYSQL_DATABASE=PipelineDatabase -e MYSQL_USER=pipeline \
-e MYSQL_PASSWORD=Run@pipelines9061 \
-d mysql/mysql-server:latest
```
## Usage
The `/workflows` directory houses all by-assay pipelines. In each subdirectory, there are several files (using RNA-seq as example):
* `rnaseq.wdl` - The main workflow file that will be run.
* `imports.zip` - A zipped directory of all imports needed for this pipeline.
* `rnaseq_inputs.json` - A JSON file with all relative input options. \

In the inputs JSON, you will find command-specific parameters (i.e. for trimming, alignment, etc.) as well as common inputs needed for every pipeline:
* `paired` - A boolean value indicating whether the experiment was paired-end (`true`) or single-end (`false`).
* `project_out_dir` - Directory where output will be (a directory is created here for every sample name).
* `fastq_dir` - Path to dorectory with all `.fastq.gz` files.
* `sampleList` - Path to your previously created `samplesheet.tsv`
* `star_index`/`BWA_index` - Path to alignment index files.
* `chromNoScaffold` - A Path to a 3-column BED file with the chromosome/contig regions you would like to keep in your output.
* `GeneAnnotationFile` - Path to `.gtf` file.
* `chromosome_sizes` - A standard BEDtools chromosome size file.
* `blacklist` - A 3-column BED file containing areas to be removed from output.

Once the configuration is complete, there are 2 options available to you to run your pipeline.
1. Run in local mode: Cromwell can run "on-the-fly", without the need to configure a server. To do this, use the `run` subcommand as such:
```
java -Dconfig.file=/path/to/cromwell_config -jar /path/to/cromwell.jar run batch_workflow.wdl -i workflow_input.json
```
2. Submit to a Cromwell server
```
==========WIP=========
```
## Scalability 
Currently, the pipeline's scalable features are implemented as `cpu` and `mem` input options in `<workflow>_inputs.json`. Each memory/compute-heavy tasks is given its own version of the variable, with `cpu` representing cores and `mem` representing memory in gigabytes. \
The aforementioned `concurrent-job-limit` variable within the Cromwell config file sets the upper limit on how many jobs, or samples, that can be processed at once.






