version 1.0

import "../../../wdl_tasks/peak_calling.wdl" as pcTasks

workflow Peakcalling {
        input {
                File sampleList
                String bam_dir
                String project_out_dir
                String? PeakcallingControl
        }
	Array[Array[String]] samples = read_tsv(sampleList)	
        scatter (sample in samples) {
		String sampleName = sample[0]
		call pcTasks.macs2 {
			input:
				bam=bam_dir+"/"+sampleName+".bam",
				sampleName=sampleName,
				control_bam=PeakcallingControl
                }
	}
}
