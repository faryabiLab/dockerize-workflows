version 1.0

# import base workflow here
import "cutAndRun.wdl" as cutandrun

#####
# The model for output file naming is as follows:
# Each sample name will be a directory, and wthin each directory will be that sample's complete collection of files
# from the entire run of the pipeline
#####

workflow batch_workflow {
	input {
		File sampleList
		String fastq_dir
		String project_out_dir
		Boolean paired
		String Aligner 
		String BWAIndex
		String PeakCaller
		String? PeakCallingControl
		Float? TopPeakFraction
		String ChromNoScaffold
		String Blacklist
		String ChromosomeSizes
		String? FastqSuffix
		Int? BamSizeFilterUpperThreshold
		Int? BamSizeFilterLowerThreshold
		Int? trimgalore_cpu = 8
		Int? trimgalore_mem = 8
		Int? bwa_cpu = 16
		Int? bwa_mem = 16
		Int? sam_to_bam_cpu = 12
		Int? sam_to_bam_mem = 8
		Int? rmScaffold_cpu = 12
		Int? rmScaffold_mem = 100
		Int? rmBlacklist_cpu = 12
		Int? rmBlacklist_mem = 25
		Int? sortbam_cpu = 12
		Int? sortBam_mem = 100
		Int? rmDuplicate_cpu = 12
		Int? rmDuplicate_mem = 25
		Int? sizeFilter_cpu = 12
		Int? sizeFilter_mem = 25
		Int? peakCalling_cpu = 12
	}

	Array[Array[String]] samples = read_tsv(sampleList)
	scatter(sample in samples) {
		call cutandrun.cut_and_run {
			input:
				trimgalore_cpu = trimgalore_cpu,
				trimgalore_mem = trimgalore_mem,
				bwa_cpu = bwa_cpu,
				bwa_mem = bwa_mem,
				sam_to_bam_cpu = sam_to_bam_cpu,
				sam_to_bam_mem = sam_to_bam_mem,
				rmScaffold_cpu = rmScaffold_cpu,
				rmScaffold_mem = rmScaffold_mem,
				rmBlacklist_cpu = rmBlacklist_cpu,
				rmBlacklist_mem = rmBlacklist_mem,
				sortbam_cpu = sortbam_cpu,
				sortBam_mem = sortBam_mem,
				rmDuplicate_cpu = rmDuplicate_cpu,
				rmDuplicate_mem = rmDuplicate_mem,
				sizeFilter_cpu = sizeFilter_cpu,
				sizeFilter_mem = sizeFilter_mem,
				peakCalling_cpu = peakCalling_cpu,
				FastqSuffix=FastqSuffix,

				sampleList=sampleList,
				fastq_dir=fastq_dir,
				sample_out_dir=project_out_dir+"/"+sample[0],
				sampleName=sample[0],
				paired=paired,
				Aligner=Aligner,
				BWAIndex=BWAIndex,
				PeakCaller=PeakCaller,
				PeakCallingControl=PeakCallingControl,
				TopPeakFraction=TopPeakFraction,
				ChromNoScaffold=ChromNoScaffold,
				Blacklist=Blacklist,
				ChromosomeSizes=ChromosomeSizes,
				size_low=BamSizeFilterLowerThreshold,
				size_high=BamSizeFilterUpperThreshold
		}	
	}
}	
