version 1.0

task count_reads_single {
	input {
		String bam
		String GeneAnnotationFile		

		String sample_name

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
                String bam
		String GeneAnnotationFile		

		String sample_name

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
