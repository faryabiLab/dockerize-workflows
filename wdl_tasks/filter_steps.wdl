version 1.0

task remove_scaffolds {
	input {
		String bam
		String chrom_no_scaff
		String sample_name
	}
	command {
		samtools view -h -L ${chrom_no_scaff} ${bam} | samtools sort - -o "${sample_name}.noScaffold.bam"
	}
	output {
		File bam_noScaffold = "${sample_name}.noScaffold.bam"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}

task remove_duplicates {
	input {
		String bam
		String sample_name

		String PicardRemoveDuplicates = "false"
		String PicardValidationStringency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
        }
	String metricsFile = "${sample_name}"+"_"+"${PicardMetricsFile}"
        command {
		java -jar $PICARD MarkDuplicates \
		M=${metricsFile} \
		O="${sample_name}.noDuplicate.bam" \
		I=${bam} \
		REMOVE_DUPLICATES=${PicardRemoveDuplicates} \
		VALIDATION_STRINGENCY=${PicardValidationStringency}
        }
        output {
 		File bam_noDuplicate = "${sample_name}.noDuplicate.bam"
        }
        runtime {
		# Picard docker image w/ samtools base
                docker: "faryabilab/picard:0.1.0"
        }
}

task remove_blacklist {
	input {
		String bam
		String blacklist
	
		String sample_name
        }
        command {
		bedtools intersect -abam ${bam} -b ${blacklist} -v > "${sample_name}.noBlacklist.bam"
        }
        output {
		File bam_noBlacklist = "${sample_name}.noBlacklist.bam"
        }
        runtime {
                docker: "faryabilab/bedtools:0.1.0"
        }
}

task sort_bam {
	input {
		String bam

		String sample_name
        }
        command {
		samtools sort ${bam} -o "${sample_name}.sorted.bam"
        }
        output {
		File bam_sorted = "${sample_name}.sorted.bam"
        }
        runtime {
                docker: "faryabilab/samtools:0.1.0"
        }

}

task index_bam {
	input {
		String bam
		
		String sample_name
	}
	command {
		samtools index -b ${bam} "${sample_name}_index.bai"
	}
	output {
		File bam_index = "${sample_name}_index.bai"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}


