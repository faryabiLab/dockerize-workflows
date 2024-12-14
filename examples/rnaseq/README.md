### RNA-seq Workflow Example
A walkthrough of setting up a workflow environment and running the workflow via Cromwell.
****
1. `rnaseq_inputs.json` modified to point to all required inputs. Inputs that are not required have their keys deleted from the input file for organization.
2. `./make_samplesheet -p -d /path/to/fastq` is run to generate `samplesheet.tsv`, with 1 sample name (`SRR26459771`). Refer to the main README for script details.
3. `options.json` modified to point to output locations for both output data and log files.
4. The workflow is run as such:
```
java -Dconfig.file=../../cromwell_configs/cromwell.docker.config -jar /path/to/cromwell.jar run -i rnaseq_inputs.json -p imports.zip -o options.json rnaseq.wdl
```
