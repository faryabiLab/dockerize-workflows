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
		String? blacklist
		String GeneAnnotationFile
		String chromosome_sizes
		# trim_galore
		Int quality = 15
		Int stringency = 5
		Float e = 0.1
		Int length = 20
		# STAR
		String outFilterType = "BySJout"
		String readFilesCommand = "zcat"
		String outSamAttributes = "Standard"
		String outFilterIntronMotifs = "RemoveNoncanonicalUnannotated"
		Int alignIntronMax = 100000
		String outSamstrandField = "intronMotif"	
		String outSAMunmapped = "Within"
		Int chimSegmentMin = 25
		Int chimJunctionOverhangMin = 25
		String outSAMtype = "BAM Unsorted"
		Int STAR_cpu=10
		Int STAR_mem=25
		# remove duplicates
		String PicardRemoveDuplicates = "true"
		String PicardValidationStringency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
	}	
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
                call rnaseq.rnaseq {
                        input:
				# Important data files for steps
                                paired=paired,
				fastq_dir=fastq_dir,
				sampleName=sample[0],
				sample_out_dir=project_out_dir+"/"+sample[0]+"/",
				star_index=star_index,
				GeneAnnotationFile=GeneAnnotationFile,
				chromosome_sizes=chromosome_sizes,
				blacklist=blacklist,
				chromNoScaffold=chromNoScaffold,
				# trim_galore 
				quality=quality,
				stringency=stringency,
				e=e,
				length=length,
				# STAR
				outFilterType=outFilterType,
				readFilesCommand=readFilesCommand,
				outSamAttributes=outSamAttributes,
				outFilterIntronMotifs=outFilterIntronMotifs,
				alignIntronMax=alignIntronMax,
				outSamstrandField=outSamstrandField,
				outSAMunmapped=outSAMunmapped,
				chimSegmentMin=chimSegmentMin,
				chimJunctionOverhangMin=chimJunctionOverhangMin,
				outSAMtype=outSAMtype,
				STAR_cpu=STAR_cpu,
				STAR_mem=STAR_mem,
				# remove duplicates
				PicardRemoveDuplicates=PicardRemoveDuplicates,
				PicardValidationStringency=PicardValidationStringency,
				PicardMetricsFile=PicardMetricsFile
                }
        }		
}
