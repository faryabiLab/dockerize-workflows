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
		String Blacklist
		String ChromosomeSizes
		String? PeakcallingControl
	}
	
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
		call chipseq.chipseq {
			input:
				paired=paired,
				fastq_dir=fastq_dir,
				sample_out_dir=project_out_dir+"/"+sample[0]+"/",
				sampleName=sample[0],
				bwa_index=BWAIndex,
				chromNoScaffold=ChromNoScaffold,
				blacklist=Blacklist,
				chromosome_sizes=ChromosomeSizes,
				PeakcallingControl=PeakcallingControl
		}
	}
}
