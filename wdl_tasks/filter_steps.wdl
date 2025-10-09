version 1.0

# testing

task sam_to_bam {
	input {
		File sam
		String sample_name
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
		Int cpu = 12
		Int mem = 8
	}
	command {
		samtools view -@ ${cpu} -h -b ${sam} > "${sample_name}.bam"
	}
	output {
		File bam = "${sample_name}.bam"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task remove_unmapped {
	input {
		#### REQUIRED
		File bam
		String sample_name
                ####
                String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
                Int cpu = 12
                Int mem = 25

		Int mem_per_thread = floor(mem / cpu)
	}
	command {
		samtools view -F 4 -f 2 -@ ${cpu} -h ${bam} | samtools sort -@ ${cpu} -m "~{mem_per_thread}G" -O bam -o "${sample_name}.noUnmapped.bam" -
	}
	output {
		File bam_noUnmapped = "${sample_name}.noUnmapped.bam"
	}
	runtime {
                docker: "${Dockerhub_Pull}"
                cpu: "${cpu}"
                mem: "${mem}"
        }
}

task remove_lowQuality {
        input {
                #### REQUIRED
                File bam
                String sample_name
                ####
                String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
                Int cpu = 12
                Int mem = 25

		Int quality_thresh = 20

                Int mem_per_thread = floor(mem / cpu)
        }
        command {
                samtools view -q ${quality_thresh} -@ ${cpu} -h ${bam} | samtools sort -@ ${cpu} -m "~{mem_per_thread}G" -O bam -o "${sample_name}.noLowQuality.bam" -
        }
        output {
                File bam_noLowQuality = "${sample_name}.noLowQuality.bam"
        }
        runtime {
                docker: "${Dockerhub_Pull}"
                cpu: "${cpu}"
                mem: "${mem}"
        }
}

task fix_mate {
        input {
                #### REQUIRED
                File bam
                String sample_name
                ####
                String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
                Int cpu = 12
                Int mem = 25

                Int mem_per_thread = floor(mem / cpu)
        }
        command {
                samtools fixmate -r -@ ${cpu} ${bam} - | samtools sort -@ ${cpu} -m "~{mem_per_thread}G" -O bam -o "${sample_name}.fix_mate.bam" -
		rm tmp.${sample_name}_filter_srt.nmsrt.bam
        }
        output {
                File bam_fixMate = "${sample_name}.fix_mate.bam"
        }
        runtime {
                docker: "${Dockerhub_Pull}"
                cpu: "${cpu}"
                mem: "${mem}"
        }
}

task remove_scaffolds {
	input {
		#### REQUIRED
		File bam
		File chrom_no_scaff
		String sample_name
		####
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
		
		Int cpu = 12
		Int mem = 25

		Int mem_per_thread = floor(mem / cpu)
	}
	command {
		samtools view -@ ${cpu} -h -L ${chrom_no_scaff} ${bam} | samtools sort -@ ${cpu} -m "~{mem_per_thread}G" -O bam -o "${sample_name}.noScaffold.bam" -
	}
	output {
		File bam_noScaffold = "${sample_name}.noScaffold.bam"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task mark_duplicates {
	input {
		#### REQUIRED
		File bam
		String sample_name
		####
		String? Dockerhub_Pull = "faryabilab/picard:0.1.0"

		String PicardValidationStringency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"

		Int cpu = 12
		Int mem = 25

        }
	String metricsFile = "${sample_name}"+"_"+"${PicardMetricsFile}"
        command {
		java -jar $PICARD MarkDuplicates \
		M=${metricsFile} \
		O="${sample_name}.dupMarked.bam" \
		I=${bam} \
		ASO=coordinate \
		REMOVE_DUPLICATES=false \
		VALIDATION_STRINGENCY=${PicardValidationStringency}
        }
        output {
 		File bam_dupMarked = "${sample_name}.dupMarked.bam"
        }
        runtime {
		# Picard docker image w/ samtools base
                docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
        }
}

task remove_duplicates_unmapped {
	input {
		File bam
		String sample_name
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
		Int cpu = 12
		Int mem = 25
	}
	command {
		samtools view -b -F 516 ${bam} > "${sample_name}.noDuplicate.bam"
	}
	output {
		File bam_noDuplicate = "${sample_name}.noDuplicate.bam"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task remove_blacklist {
	input {
		#### REQUIRED
		File bam
		File? blacklist
		String sample_name	
		####
		String? Dockerhub_Pull = "faryabilab/bedtools:0.1.0"

		Int cpu = 12
		Int mem = 25
        }
        command {
		bedtools intersect -abam ${bam} -b ${blacklist} -v > "${sample_name}.noBlacklist.bam"	
        }
        output {
		File bam_noBlacklist = "${sample_name}.noBlacklist.bam"
        }
        runtime {
                docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
        }
}

task sort_bam {
	input {
		File bam
		String sample_name
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"

		Int cpu = 12
		Int mem = 100

		Int mem_per_thread = floor(mem / cpu)
        }
        command {
		samtools sort -@ ${cpu} -m "${mem_per_thread}G" ${bam} -o "${sample_name}.sorted.bam"
        }
        output {
		File bam_sorted = "${sample_name}.sorted.bam"
        }
        runtime {
                docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
        }
}

task sort_bam_name {
        input {
                File bam
                String sample_name
                String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"

                Int cpu = 12
                Int mem = 100

                Int mem_per_thread = floor(mem / cpu)
        }
        command {
                samtools sort -n -@ ${cpu} -m "${mem_per_thread}G" ${bam} -o "${sample_name}.sortedByName.bam"
        }
        output {
                File bam_sorted = "${sample_name}.sortedByName.bam"
        }
        runtime {
                docker: "${Dockerhub_Pull}"
                cpu: "${cpu}"
                mem: "${mem}"
        }
}

task index_bam {
	input {
		File bam
		String sample_name
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
		
		Int cpu = 12
		Int mem = 25
	}
	command {
		samtools index -@ ${cpu} -b ${bam} "${sample_name}_index.bai"
	}
	output {
		File bam_index = "${sample_name}_index.bai"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task size_filter_bam {
	input {
		File bam
		String sample_name
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
		
		Int? threshold_low
		Int? threshold_hi
		
		Int cpu = 12
		Int mem = 5
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
		File? low = "${sample_name}.lowThreshold_GreaterThan${threshold_low}bp.bam"
		File? hi = "${sample_name}.hiThreshold_LessThan${threshold_hi}bp.bam"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task filter_discordant_pairs {
	input {
		File bam
		String sample_name
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"
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
		File bam_pairedReads = "${sample_name}.pairedReads.bam"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
