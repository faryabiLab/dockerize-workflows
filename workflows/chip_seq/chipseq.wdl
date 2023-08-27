version 1.0

import "../../wdl_tasks/trim_galore.wdl" as trimTasks
import "../../wdl_tasks/bwa.wdl" as bwaTasks
import "../../wdl_tasks/filter_steps.wdl" as filterTasks
import "../../wdl_tasks/feature_count.wdl" as quantTasks
import "../../wdl_tasks/make_bigWig.wdl" as bwtasks

workflow chipseq {
	input {
		Boolean paired
                String fastq_dir
		String sample_out_dir
                String sampleName       
                String bwa_index
                String chromNoScaffold
		String blacklist
		String chromosome_sizes
	}

	call trimTasks.fastqc_trim {
		input:
			fastq_dir=fastq_dir,
			sample_out_dir=sample_out_dir,
			sampleName=sampleName,
			paired=paired
	}
}
