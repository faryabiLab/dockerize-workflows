version 1.0

task count_reads_single {
	input {
		File bam
		File GeneAnnotationFile		

		String out_dir
		String sample_name

		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"
	}
	out = "${out_dir}"+"/"+"${sample_name}"+".counts"
	command {
		featureCounts \	
			-t ${AttributeType} \
			-g ${GTFAttributeType} \
			-s ${Stranded} \
			-a ${GeneAnnotationFile} \
			-o ${out} \
			${bam}
	}
	output {
		File gene_counts = ${out}
	}
	runtime {
		docker: "faryabilab/subread:0.1.0"
	}
}

task count_reads_paired {
        input {
                File bam
		File GeneAnnotationFile		

		String out_dir
		String sample_name

		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"
        }
        out = "${out_dir}"+"/"+"${sample_name}"+".counts"
        command {
		featureCounts \
			-p \
                        -t ${AttributeType} \
                        -g ${GTFAttributeType} \
                        -s ${Stranded} \
                        -a ${GeneAnnotationFile} \
                        -o ${out} \
                        ${bam}
        }
        output {
		File gene_counts = ${out}
        }
        runtime {
                docker: "faryabilab/subread:0.1.0"
        }
}
