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

                Int qualityCutoff = 15
                Int stringencyCutoff = 5
                Float errorRate = 0.1
                Int lengthCutoff = 20
                # STAR args             
                String star_index
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
                # Remove scaffold args
                String chromNoScaffold
		# Remove duplicates args
		String PicardRemoveDuplicates = "false"
		String PicardValidationStringency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
		# Remove blacklist args
		File blacklist
		# Sort bam args
		# Index bam args
		# Quantify args
		String GeneAnnotationFile
		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"
		# Make BigWig args
		String chromosome_sizes
	}
	
	Array[Array[String]] samples = read_tsv(sampleList)
	if (paired) {
		scatter (sample in samples) {
                	call rnaseq.rnaseq {
                        	input:
                                	paired=paired,
					fastq_dir=fastq_dir,
					sampleName=sample[0],
					sample_out_dir=project_out_dir+"/"+sample[0]+"/",
					qualityCutoff=qualityCutoff,
					stringencyCutoff=stringencyCutoff,
                			errorRate=errorRate,
                			lengthCutoff=lengthCutoff,
					star_index=star_index,
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
					chromNoScaffold=chromNoScaffold,
					PicardRemoveDuplicates=PicardRemoveDuplicates,
					PicardValidationStringency=PicardValidationStringency,
					PicardMetricsFile=PicardMetricsFile,
					GeneAnnotationFile=GeneAnnotationFile,
					AttributeType=AttributeType,
					GTFAttributeType=GTFAttributeType,
					Stranded=Stranded,
					chromosome_sizes=chromosome_sizes,
					blacklist=blacklist

                	}
        	}	
	}
	
		
}
	







