version 1.0

task remove_scaffolds {
	input {
		File bam
		File chrom_no_scaff
		String out_dir
		String sample_name
	}
	String out = "${out_dir}"+"/"+"${sample_name}"+".noScaffold.bam"
	command {
		samtools view -h -L ${chrom_no_scaff} ${bam} | samtools sort - -o ${out}
	}
	output {
		File j = ${out}
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}

task remove_duplicates {
	input {
		File bam
		String out_dir
		String sample_name

		String PicardRemoveDuplicates = "false"
		String PicardValidationStrignency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
        }
	out = "${out_dir}"+"/"+"${sample_name}"+".noDuplicate.bam"
	metricsFile = "${sample_name}"+"_"+"${PicardMetricsFile}"
        command {
		picard MarkDuplicates \
		M=${metricsFile} \
		O=${out} \
		I=${bam_noScaffold} \
		REMOVE_DUPLICATES=${PicardRemoveDuplicates} \
		VALIDATION_STRINGENCY=${PicardValidationStringency}
        }
        output {
 		File bam_noscaff = ${out}
        }
        runtime {
		# Picard docker image w/ samtools base
                docker: "faryabilab/picard:0.1.0"
        }
}

task remove_blacklist {
	input {
		File bam
		File blacklist
	
		String out_dir
		String sample_name
        }
	out = "${out_dir}"+"/"+"${sample_name}"+".noBlacklist.bam"
        command {
		bedtools intersect -abam ${bam} -b ${blacklist} -v > ${out}
        }
        output {
		File bam_noBlacklist = ${out}
        }
        runtime {
                docker: "faryabilab/bedtools:0.1.0"
        }
}

task sort_bam {
	input {
		File bam

		String out_dir
		String sample_name
        }
	out = ${out_dir}+"/"+"${sample_name}"+".sorted.bam"
        command {
		samtools sort ${bam} -o ${out}
        }
        output {
		File bam_sorted = ${bam}
        }
        runtime {
                docker: "faryabilab/samtools:0.1.0"
        }

}

task index_bam {
	input {
		File bam

		String out_dir
		String sample_name
	}
	out = "${out+dir}"+"/"+"${sample_name}"+"_index.bai"
	command {
		samtools index -b ${bam} ${out}
	}
	output {
		File bam_index = ${out}
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}


