version 1.0

import "../../../wdl_tasks/trim_galore.wdl" as trimTasks

workflow RNAseq {
	input {
		String sampleList
		String fastq_dir
		Boolean paired       
		String? fastq_suffix
	}	
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
		String sampleName=sample[0]
		call trimTasks.fastqc_trim {
			input:
				fastq_dir=fastq_dir,
				sampleName=sampleName,
				paired=paired,
				fastq_suffix=fastq_suffix
		}
			
        }
	output {
		Array[File?] fastqc = fastqc_trim.out_html
	}
}
