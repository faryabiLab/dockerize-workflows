version 1.0

task count_reads_single {
	input {
		#### REQUIRED
		File bam
		String GeneAnnotationFile		
		String sample_name
		####

		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"

		Int cpu = 12
		Int mem = 25
	}
	command {
		featureCounts \
			-T ${cpu} \
			-t ${AttributeType} \
			-g ${GTFAttributeType} \
			-s ${Stranded} \
			-a ${GeneAnnotationFile} \
			-o "${sample_name}.counts" \
			${bam}
	}
	output {
		File gene_counts = "${sample_name}.counts"
	}
	runtime {
		docker: "faryabilab/subread:0.1.0"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task count_reads_paired {
        input {
		#### REQUIRED
                File bam
		String GeneAnnotationFile		
		String sample_name
		####

		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"

		Int cpu = 12
		Int mem = 25
        }
        command {
		featureCounts \
		-T ${cpu} \
		-p \
		-t ${AttributeType} \
		-g ${GTFAttributeType} \
		-s ${Stranded} \
		-a ${GeneAnnotationFile} \
		-o "${sample_name}.counts" \
		${bam}
        }
        output {
		File gene_counts = "${sample_name}.counts"
        }
        runtime {
                docker: "faryabilab/subread:0.1.0"
		cpu: "${cpu}"
		mem: "${mem}"
        }
}

task quantifyCoverage {
	input {
		#### REQUIRED
		File bam
		String peaks
		String chromSizes
		String sample_name
		####

		Int cpu = 12
		Int mem = 25
	}
	command {
		bedtools coverage \
		-a ${peaks} -b ${bam} \
		-g ${chromSizes} \
		> "${sample_name}.coverage.counts"
	}
	output {
		File out_counts = "${sample_name}.coverage.counts"
	}
	runtime {
		docker: "faryabilab/bedtools:0.1.0"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
