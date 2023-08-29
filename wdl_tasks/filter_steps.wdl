version 1.0

task sam_to_bam {
	input {
		String sam
		String sample_name
	}
	command {
		samtools view -h -b ${sam} > "${sample_name}.bam"
	}
	output {
		String bam = "${sample_name}.bam"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}

task remove_scaffolds {
	input {
		#### REQUIRED
		String bam
		String chrom_no_scaff
		String sample_name
		####
	}
	command {
		samtools view -h -L ${chrom_no_scaff} ${bam} | samtools sort - -o "${sample_name}.noScaffold.bam"
	}
	output {
		String bam_noScaffold = "${sample_name}.noScaffold.bam"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}

task remove_duplicates {
	input {
		#### REQUIRED
		String bam
		String sample_name
		####

		String PicardRemoveDuplicates = "true"
		String PicardValidationStringency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
        }
	String metricsFile = "${sample_name}"+"_"+"${PicardMetricsFile}"
        command {
		java -jar $PICARD MarkDuplicates \
		M=${metricsFile} \
		O="${sample_name}.noDuplicate.bam" \
		I=${bam} \
		ASO=coordinate \
		REMOVE_DUPLICATES=${PicardRemoveDuplicates} \
		VALIDATION_STRINGENCY=${PicardValidationStringency}
        }
        output {
 		String bam_noDuplicate = "${sample_name}.noDuplicate.bam"
        }
        runtime {
		# Picard docker image w/ samtools base
                docker: "faryabilab/picard:0.1.0"
        }
}

task remove_blacklist {
	input {
		#### REQUIRED
		String bam
		String? bam2
		String blacklist
		String sample_name
		Boolean? aligner
		####
        }
        command {
		if [ "${aligner}" == "bowtie2" ]; then
			bedtools intersect -abam ${bam} -b ${blacklist} -v > "${sample_name}.noBlacklist.bam"
		elif [ "${aligner}" == "bwa" ]; then
			bedtools intersect -abam ${bam2} -b ${blacklist} -v > "${sample_name}.noBlacklist.bam"
			
		else
			bedtools intersect -abam ${bam} -b ${blacklist} -v > "${sample_name}.noBlacklist.bam"
		fi
        }
        output {
		String bam_noBlacklist = "${sample_name}.noBlacklist.bam"
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
		String bam_sorted = "${sample_name}.sorted.bam"
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
		String bam_index = "${sample_name}_index.bai"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}

task size_filter_bam {
	input {
		String bam
		String sample_name
		
		Int? threshold_low = 120
		Int? threshold_hi = 150
	}
	command {
		if [ ! -z "${threshold_low}" ] && [ ! -z "${threshold_hi}" ]; then
			samtools view \
			-e 'length(seq)<=${threshold_hi}' 
			-O BAM \
			-o "${sample_name}.hiThreshold_LessThan${threshold_hi}bp.bam" \
			${bam}
			samtools view \
			-e 'length(seq)>${threshold_low}' \
			-O BAM \
			-o "${sample_name}.lowThreshold_GreaterThan${threshold_low}bp.bam" \
			${bam}
		elif [ ! -z "${threshold_low}" ] && [ -z "${threshold_hi}" ]; then
			samtools view \
			-e 'length(seq)>=${threshold_low}' \
			-O BAM \
			-o "${sample_name}.lowThreshold_GreaterThan${threshold_low}bp.bam" \
			${bam}
		elif [ -z "${threshold_low}" ] && [ ! -z "${threshold_hi}" ]; then
			samtools view \
			-e 'length(seq)<=${threshold_hi}' \
			-O BAM \
			-o "${sample_name}.hiThreshold_LessThan${threshold_hi}bp.bam" \
			${bam}
		else
			echo ""
		fi
	}
	output {
		File? low = "${sample_name}.lowThreshold_${threshold_low}bp.bam"
		File? hi = "${sample_name}.hiThreshold_${threshold_hi}bp.bam"
	}
	runtime {
		docker: 'faryabilab/samtools:0.1.0'
	}
}
