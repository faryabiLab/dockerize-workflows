version 1.0

task count_reads_single {
	input {
		#### REQUIRED
		String bam
		String GeneAnnotationFile		
		String sample_name
		####

		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"
	}
	command {
		featureCounts \	
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
	}
}

task count_reads_paired {
        input {
		#### REQUIRED
                String bam
		String GeneAnnotationFile		
		String sample_name
		####

		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"
        }
        command {
		featureCounts \
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
        }
}

task quantifyCoverage {
	input {
		#### REQUIRED
		String bam
		String peaks
		String chromSizes
		String sample_name
		####
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
	}
}












