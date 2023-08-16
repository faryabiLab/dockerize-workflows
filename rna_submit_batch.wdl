version 1.0

import "rnaseq.wdl" as rnaseq

workflow batch_workflow {
	input {
		File sampleList
		String fastq_dir
		String project_out_dir = "."

		Boolean paired

		String trim_out_dir = "01.trim_fastqc"
                Int qualityCutoff = 15
                Int stringencyCutoff = 5
                Float errorRate = 0.1
                Int lengthCutoff = 20
                # STAR args             
                String star_index
                String star_out_dir = "02.alignment"
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
                File chromNoScaffold
		String removeScaffold_out_dir = "03.removeScaffold"
		# Remove duplicates args
		String removeDuplicate_out_dir = "04.removeDuplicate"
		String PicardRemoveDuplicates = "false"
		String PicardValidationStringency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
		# Remove blacklist args
		String removeBlacklist_out_dir = "05.removeBlacklist"
		File blacklist
		# Sort bam args
		String sortBam_out_dir = "06.sortedIndexedBam"
		# Index bam args
		String indexBam_out_dir = "06.sortedIndexedBam"
		# Quantify args
		String counts_out_dir = "07.counts"
		File GeneAnnotationFile
		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"
		# Make BigWig args
		String bw_out_dir = "08.bigWigs"
		File chromosome_sizes
	}
	
	Array[Array[String]] samples = read_tsv(sampleList)
	if (paired) {
		scatter (sample in samples) {
                	call rnaseq.rnaseq {
                        	input:
                                	paired=paired,
					fastq1="${fastq_dir}"+"/"+sample[0]+"_R1.fastq.gz",
					fastq2="${fastq_dir}"+"/"+sample[0]+"_R2.fastq.gz",
					sampleName=sample[0],
					trim_out_dir="${project_out_dir}/${trim_out_dir}",
					qualityCutoff=qualityCutoff,
					stringencyCutoff=stringencyCutoff,
                			errorRate=errorRate,
                			lengthCutoff=lengthCutoff,
					star_index=star_index,
                			star_out_dir="${project_out_dir}/${star_out_dir}",
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
					removeScaffold_out_dir="${project_out_dir}/${removeScaffold_out_dir}",
					removeDuplicate_out_dir="${project_out_dir}/${removeDuplicate_out_dir}",
					PicardRemoveDuplicates=PicardRemoveDuplicates,
					PicardValidationStringency=PicardValidationStringency,
					PicardMetricsFile=PicardMetricsFile,
					removeBlacklist_out_dir="${project_out_dir}/${removeBlacklist_out_dir}",
					sortBam_out_dir="${project_out_dir}/${sortBam_out_dir}",
					indexBam_out_dir="${project_out_dir}/${indexBam_out_dir}",
					counts_out_dir="${project_out_dir}/${counts_out_dir}",
					GeneAnnotationFile=GeneAnnotationFile,
					AttributeType=AttributeType,
					GTFAttributeType=GTFAttributeType,
					Stranded=Stranded,
					bw_out_dir="${project_out_dir}/${bw_out_dir}",
					chromosome_sizes=chromosome_sizes,
					blacklist=blacklist

                	}
        	}	
	}
	
		
}
	







