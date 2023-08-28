version 1.0

import "../../wdl_tasks/trim_galore.wdl" as trimTasks
import "../../wdl_tasks/bowtie2.wdl" as bowtieTask
import "../../wdl_tasks/filter_steps" as filterTasks
import "../../wdl_tasks/feature_count" as quantTasks
import "../../wdl_tasks/peak_calling.wdl" as pcTasks
import "../../wdl_tasks/make_bigWig.wdl" as bwTasks

workflow cut_and_run {
	input {
		File sampleList
		String fastq_dir
		String sample_out_dir
		String sampleName
		Boolean paired
		String Aligner 
		String? BWAIndex
		String? BowtieIndex
		String? PeakCaller
		String PeakCallingControl
		String ChromNoScaffold
		String Blacklist
		String ChromosomeSizes
	}
	# call trim_galore
	call trimTasks.fastqc_trim {
		input:
			fastq_dir=fastq_dir,
			
			
	}
}
