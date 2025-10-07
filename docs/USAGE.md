## Sample Sheet
Every workflow run is accomapnied by a `samplesheet.tsv`. It is a 3-column tab-separated file containing a **sample name**, the **path to the associated fastq file**, and the **paired read file (i.e. R2)** if applicable.  \

To make a sample sheet, use `utils/make_samplesheet.sh`:
* Paired-end (named as `*_<R1/R2>.fastq.gz`): `./make_samplesheet.sh -p -s R -d /path/to/fastq_dir`
* Paired-end (named as `*_<1/2>.fastq.gz`): `./make_samplesheet.sh -p -s N -d /path/to/fastq_dir`
* Single-end: `./make_samplesheet.sh -d /path/to/fastq_dir` 

This will create a `samplesheet.tsv` in your current working directory.

## Reference Genome
The workflow input requires that the directory holding reference genome index files be zipped. This can be done with the script:
```
utils/zip_reference.sh -d <dir>
```

## Call-Caching
A desirable feature for many workflows is the ability to identify steps, or "jobs", that have already been run and essentially "skip" that step, saving both time and computing resources. In this reporisoty, Cromwell is configured to use a local MySQL database (which runs from a Docker container) called `PipelineDatabase`, to store intermediate files to determine the steps that have been run. If you want to use a different SQL, or any, database, the `database` block must be edited accordingly in the Cromwell configuration file, else, the config will attempt to log in to the MySQL server specified. By default, the `database` block is configured to work with a MySQL Docker container that can be initialized by running the script:
```
./workflows/start_database.sh
```
This will create `index.tar.gz` in the current working directory. Use this folder when passing directories as parameters (below).

## Options file
Cromwell allows the user to pass in an `options` file via the `-o/--options` argument, which grants control over the pipeline's output locations. This file's format is that of standard JSON, and the keys made readily available in this repository are:
* `final_workflow_outputs_dir` - Directory to which the final workflow outputs will be copied
* `final_call_logs_dir` - Directory to which the final workflow logs will be copied
There are more options available, you can find them in the Cromwell documentation.

## Docker Containers
By default, the pipelines are configured to run each step within a specified docker container, all of which are hosted on the [Faryabi Lab Dockerhub](https://hub.docker.com/). Workflow input files contain a `Dockerhub_Pull` variable for each step that can be used to point to another Docker container, although this could be risky.

## Cromwell Configuration
Within the `cromwell_configs` directory is a Cromwell config file which instructs the engine to use Docker as backend. The only option that can be tweaked in this file is the `concurrent_job_limit` variable: this controls the number of jobs that can be running at once (default = 10).

## Running a Workflow
### Files
The `/workflows` directory houses all by-assay pipelines. In each subdirectory, there are several files:
* `*.wdl` - The main workflow file that will be run.
* `imports.zip` - A zipped directory of all task imports needed for this pipeline. This is important for running your pipeline successfully. 
* `*_inputs.json` - A JSON file with all relative input options.
* `options.json` - A JSON file with paths to output locations for workflow results and logs. 

### Input Parameters
In the inputs JSON, you will find command-specific parameters (i.e. for trimming, alignment, etc.) as well as common inputs needed for every pipeline:
* `paired` - A boolean value indicating whether the experiment was paired-end (`true`) or single-end (`false`).
* `project_out_dir` - Directory where output will be (a directory is created here for every sample name).
* `sampleList` - Path to your previously created `samplesheet.tsv`
* `star_index`/`BWA_index` - Path to the zipped index directory.
* `chromNoScaffold` - A Path to a 3-column BED file with the chromosome/contig regions you would like to keep in your output.
* `GeneAnnotationFile` - Path to `.gtf` file.
* `chromosome_sizes` - A standard BEDtools chromosome size file.
* `blacklist` (Optional) - A 3-column BED file containing areas to be removed from output.
### Executing
Once the configuration is complete, there are 2 options available to you to run your pipeline: 

**1. Run in local mode**
* Cromwell can run "on-the-fly", without the need to configure a server. To do this, use `workflows/run_workflow.sh`
```
./run_workflow \
  -j <path to cromwell .jar> \
  -c <cromwell config file> \
  -i <input JSON> \
  -o <cromwell options JSON> \
  -p <zipped imports file> \
  -w <.wdl workflow file>
```
* This will start a Cromwell instance, run your workflow, and then exit upon completion.

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
