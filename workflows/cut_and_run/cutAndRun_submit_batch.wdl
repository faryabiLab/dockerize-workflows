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
		String? BWAIndex
		String? BowtieIndex
		String PeakCaller
		String? PeakCallingControl
		String ChromNoScaffold
		String Blacklist
		String ChromosomeSizes
		Int? BamSizeFilterUpperThreshold = 120
		Int? BamSizeFilterLowerThreshold = 150
	}

	Array[Array[String]] samples = read_tsv(sampleList)
	scatter(sample in samples) {
		call cutandrun.cut_and_run {
			input:
				sampleList=sampleList,
				fastq_dir=fastq_dir,
				sample_out_dir=project_out_dir+"/"+sample[0],
				sampleName=sample[0],
				paired=paired,
				Aligner=Aligner,
				BWAIndex=BWAIndex,
				BowtieIndex=BowtieIndex,
				PeakCaller=PeakCaller,
				PeakCallingControl=PeakCallingControl,
				ChromNoScaffold=ChromNoScaffold,
				Blacklist=Blacklist,
				ChromosomeSizes=ChromosomeSizes,
				size_low=BamSizeFilterLowerThreshold,
				size_high=BamSizeFilterUpperThreshold
		}	
	}
}	

















