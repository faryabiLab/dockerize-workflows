version 1.0

task remove_scaffolds {
	input {
		File bam
		File chrom_no_scaff
		String out_dir = "02.alignment"
	}
	prefix = basename(bam, ".bam")
	out = "${out_dir}"+"/"+"${prefix}"+".noScaffold.bam"
	command {
		samtools view -h -L ${chrom_no_scaff} ${raw_bam} | samtools sort - -o ${out}
	}
	output {
		File bam_sorted = ${out}
	}
	runtime {
		docker: "faryabilab/samtools:0.1.0"
	}
}

task remove_duplicates {
	input {
		File bam
		String out_dir = "02.alignment"

		String PicardRemoveDuplicates = "false"
		String Picard CalidationStrignency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
        }
	prefix = basename(bam_noScaffold, ".bam")
	out = "${out_dir}"+"/"+"${out_prefix}"+".noDuplicate.bam"
	metricsFile = prefix+"_"+"${PicardMetricsFile}"
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
	
		String out_dir = "02.alignment"
        }
	prefix = basename(bam, ".bam")
	out = "${out_dir}"+"/"+"${prefix}+".noBlacklist.bam"
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
	input{
		File bam

		String out_dir
        }
	prefix = basename(bam, ".bam")
	out = ${out_dir}+"/"+"${prefix}"+".sorted.bam"
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
	}
	prefix = basename(bam, ".bam")
	out = "${out+dir}"+"/"+"${prefix}"+"_index.bai"
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


