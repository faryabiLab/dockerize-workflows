# Dockerized Genomic Processing Pipelines
Genomic data processing pipelines written in Workflow Definition Language (WDL), making use of Docker containers to ensure reproducability. The Docker containers used have been curated for usage in the Faryabi Lab processing pipelines, and can be found on [Dockerhub](https://hub.docker.com/u/faryabilab).

## Setup
Clone this repository to your local machine
```
git clone https://github.com/faryabiLab/dockerize-workflows
```
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
#### Options file
Cromwell allows the user to pass in an `options` file via the `-o/--options` argument, which grants control over the pipeline's output locations. This file's format is that of standard JSON, and the keys made readily available in this repository are:
* `final_workflow_outputs_dir` - Directory to which the final workflow outputs will be copied
* `final_call_logs_dir` - Directory to which the final workflow logs will be copied
There are more options available, you can find them in the Cromwell documentation.
#### Docker Containers
By default, the pipelines are configured to run each step within a specified docker container, all of which are hosted on the [Faryabi Lab Dockerhub](https://hub.docker.com/). To learn how to use different docker containers for running workflow steps, refer to the [Official Cromwell Documentation](https://cromwell.readthedocs.io/en/stable/).
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
### Inputs
The `/workflows` directory houses all by-assay pipelines. In each subdirectory, there are several files (using RNA-seq as example):
* `rnaseq.wdl` - The main workflow file that will be run.
* `imports.zip` - A zipped directory of all imports needed for this pipeline. This is important for running your pipeline successfully. 
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
### Running a Workflow
Once the configuration is complete, there are 2 options available to you to run your pipeline: 

**1. Run in local mode**
* Cromwell can run "on-the-fly", without the need to configure a server. To do this, use the `run` subcommand as such:
```
java \
  -Dconfig.file=/path/to/cromwell_config \
  -jar /path/to/cromwell.jar run \
  -i workflow_input.json \
  -o options.json \
  -p imports.zip \
  workflow.wdl
```
* This will startg a Cromwell instance, run your workflow, and then exit upon completion.

**2. Submit to a Cromwell server**
* First, if a Cromwell server isn't already running, start one:
```
java -Dconfig.file=/path/to/cromwell_config -jar /path/to/cromwell.jar server
```
* Then, submit your workflow to the server via the `submit` command along with the host address, `imports.zip`, `inputs.json`, and your workflow's WDL file:
```
java \
  -Dconfig.file=/path/to/cromwell_config \
  -jar /path/to/cromwell.jar submit \
  -h http://<host ip>:<port> \
  -p imports.zip \
  -o options.json \
  -i inputs.json \
   workflow.wdl
```
* Navigate to `0.0.0.0:8000` to access the GUI, where different aspects of the workflow can be monitored via the `REST` API.
  * A useful tool is the timing diagram tool which displays a graphical representation of the workflow's progress, broken down by sample and job. It can be accessed with the `REST` API via `0.0.0.0:5200/api/workflows/v1/<workflow ID>/timing` \

Of course, you can find out more about Cromwell's capabilities in the [Official Cromwell Documentation](https://cromwell.readthedocs.io/en/stable/).


