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
		
	}
}
