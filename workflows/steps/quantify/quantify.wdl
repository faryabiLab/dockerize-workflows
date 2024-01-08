version 1.0

import "../../../wdl_tasks/feature_count.wdl" as quantTasks

workflow Quantify {
        input {
                File sampleList
		File regions
		String ChromosomeSizes
                String bam_dir
                String project_out_dir
        }
        Array[Array[String]] samples = read_tsv(sampleList)
        scatter (sample in samples) {
                String sampleName = sample[0]
                call quantTasks.quantifyCoverage {
                        input:
				peaks=regions,
                                bam=bam_dir+"/"+sampleName+".bam",
                                sample_name=sampleName,
				chromSizes=ChromosomeSizes
                }
        }
} 
