version 1.0

import "chipseq.wdl" as chipseq

#####
# The model for output file naming is as follows:
# Each sample name will be a directory, and wthin each directory will be that sample's complete collection of files
# from the entire run of the pipeline
#####

workflow batch_workflow {
	input {
		File sampleList
		String fastq_dir
		String project_out_dir = "."
		Boolean paired
		String BWAIndex
		String ChromNoScaffold
		Stirng Blacklist
		String ChromosomeSizes
	}
	
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
		call chipseq.chipseq {
			input:
		}
	}
}
