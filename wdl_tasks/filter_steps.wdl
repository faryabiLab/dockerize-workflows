version 1.0

task sam_to_bam {
	input {
		String sam
		String sample_name

		Int cpu = 12
		Int mem = 8
	}
	command {
		samtools view -@ ${cpu} -h -b ${sam} > "${sample_name}.bam"
	}
	output {
		String bam = "${sample_name}.bam"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
		cpu: "${cpu}"
		#memory: "${mem}"
	}
}

task remove_scaffolds {
	input {
		#### REQUIRED
		String? bam
		String? bam2
		String chrom_no_scaff
		String sample_name
		####
		
		Int cpu = 12
		Int mem = 25
	}
	command {
		if [ ! -z "${bam}" ]; then 
			samtools view -@ ${cpu} -h -L ${chrom_no_scaff} ${bam} | samtools sort -@ ${cpu} -m "${mem}G" -O bam -o "${sample_name}.noScaffold.bam" -
		else
			samtools view -@ ${cpu} -h -L ${chrom_no_scaff} ${bam2} | samtools sort -@ ${cpu} -m "${mem}G" -O bam -o "${sample_name}.noScaffold.bam" -
		fi
	}
	output {
		String bam_noScaffold = "${sample_name}.noScaffold.bam"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
		cpu: "${cpu}"
		#memory: "${mem}"
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

		Int cpu = 12
		Int mem = 25
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
		cpu: "${cpu}"
		#memory: "${mem}"
        }
}

task remove_blacklist {
	input {
		#### REQUIRED
		String bam
		String? bam2
		String blacklist
		String sample_name	
		####

		Int cpu = 12
		Int mem = 25
        }
        command {
		bedtools intersect -abam ${bam} -b ${blacklist} -v > "${sample_name}.noBlacklist.bam"	
        }
        output {
		String bam_noBlacklist = "${sample_name}.noBlacklist.bam"
        }
        runtime {
                docker: "faryabilab/bedtools:0.1.0"
		cpu: "${cpu}"
		#memory: "${mem}"
        }
}

task sort_bam {
	input {
		String bam
		String sample_name

		Int cpu = 12
		Int mem = 100
        }
        command {
		samtools sort -@ ${cpu} -m "${mem}G" ${bam} -o "${sample_name}.sorted.bam"
        }
        output {
		String bam_sorted = "${sample_name}.sorted.bam"
        }
        runtime {
                docker: "faryabilab/samtools:0.1.0"
		cpu: "${cpu}"
		#memory: "${mem}"
        }
}

task index_bam {
	input {
		String bam
		String sample_name
		
		Int cpu = 12
		Int mem = 25
	}
	command {
		samtools index -@ ${cpu} -b ${bam} "${sample_name}_index.bai"
	}
	output {
		String bam_index = "${sample_name}_index.bai"
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
		cpu: "{cpu}"
		#memory: "${mem}"
	}
}

task size_filter_bam {
	input {
		String bam
		String sample_name
		
		Int? threshold_low
		Int? threshold_hi
		
		Int cpu = 12
		Int mem = 25
	}
	command {
		if [ ! -z "${threshold_low}" ] && [ ! -z "${threshold_hi}" ]; then
			samtools view \
			-e 'tlen <= ${threshold_hi} && tlen >= -${threshold_hi}' \
			-@ ${cpu} \
			-O BAM \
			-o "${sample_name}.hiThreshold_LessThan${threshold_hi}bp.bam" \
			${bam}
			samtools view \
			-e 'tlen >= ${threshold_low} && tlen <= -${threshold_low}' \
			-@ ${cpu} \
			-O BAM \
			-o "${sample_name}.lowThreshold_GreaterThan${threshold_low}bp.bam" \
			${bam}
		elif [ ! -z "${threshold_low}" ] && [ -z "${threshold_hi}" ]; then
			samtools view \
			-e 'tlen >= ${threshold_low} && tlen <= -${threshold_low}' \
			-@ ${cpu} \
			-O BAM \
			-o "${sample_name}.lowThreshold_GreaterThan${threshold_low}bp.bam" \
			${bam}
		elif [ -z "${threshold_low}" ] && [ ! -z "${threshold_hi}" ]; then
			samtools view \
			-e 'tlen <= ${threshold_hi} && tlen >= -${threshold_hi}' \
			-@ ${cpu} \
			-O BAM \
			-o "${sample_name}.hiThreshold_LessThan${threshold_hi}bp.bam" \
			${bam}
		else
			echo ""
		fi
	}
	output {
		String? low = "${sample_name}.lowThreshold_GreaterThan${threshold_low}bp.bam"
		String? hi = "${sample_name}.hiThreshold_LessThan${threshold_hi}bp.bam"
	}
	runtime {
		docker: 'faryabilab/samtools:0.1.0'
		cpu: "${cpu}"
		#memory: "${mem}"
	}
}

task filter_discordant_pairs {
	input {
		String bam
		String sample_name
		Int cpu = 12
		Int mem = 25
	}
	command {
		samtools view \
		-@ ${cpu} \
		-O BAM \
		-f 2 \
		-o "${sample_name}.pairedReads.bam" \
		${bam}
	}
	output {
		String bam_pairedReads = "${sample_name}.pairedReads.bam"
	}
	runtime {
		docker: 'faryabilab/samtools:0.1.0'
		cpu: "${cpu}"
	}
}
