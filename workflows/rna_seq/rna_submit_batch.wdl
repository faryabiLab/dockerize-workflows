version 1.0

import "rnaseq.wdl" as rnaseq

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
                String star_index
                String chromNoScaffold
		String blacklist
		String GeneAnnotationFile
		String chromosome_sizes
	}
	
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
                call rnaseq.rnaseq {
                        input:
                                paired=paired,
				fastq_dir=fastq_dir,
				sampleName=sample[0],
				sample_out_dir=project_out_dir+"/"+sample[0]+"/",
				star_index=star_index,
				GeneAnnotationFile=GeneAnnotationFile,
				chromosome_sizes=chromosome_sizes,
				blacklist=blacklist,
				chromNoScaffold=chromNoScaffold
                }
        }	
		
}
	







