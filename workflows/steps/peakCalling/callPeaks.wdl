version 1.0

import "../../../wdl_tasks/peak_calling.wdl" as pcTasks

workflow Peakcalling {
        input {
                File sampleList
                String bam_dir
                String project_out_dir
                String? PeakcallingControl
		String? peakCaller="macs2"
		String? top_peak_fraction
        }
	Array[Array[String]] samples = read_tsv(sampleList)	
        scatter (sample in samples) {
		String sampleName = sample[0]
		if ("${peakCaller}" == "macs2") {	
			call pcTasks.macs2 {
				input:
					bam=bam_dir+"/"+sampleName+".bam",
					sampleName=sampleName,
					control_bam=PeakcallingControl
                	}
		}
		if ("${peakCaller}" == "seacr") {
			call pcTasks.bamToBedgraph {
				input: 
					bam=bam_dir+"/"+sampleName+".bam",
					sampleName=sampleName,
					type="all"
			}
			call pcTasks.SEACR {
				input:
					sampleName=sampleName,
					bedgraph = bamToBedgraph.seacr_bg,
					control_bedgraph=PeakcallingControl,
					top_peak_fraction=top_peak_fraction,
					type="all"
			}
		}
	}
}
