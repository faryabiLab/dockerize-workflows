version 1.0

task count_reads_single {
	input {
		File bam
		File GeneAnnotationFile		

		String out_dir = "03.counts"
		String AttributeType = "exon"
		String GTFAttributeType = "gene_id"
		String Stranded = "1"
	}
	prefix = basename(bam, ".bam")
	out = "${out_dir}"+"/"+"${prefix}"+".counts"
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

                String out_dir = "03.counts"
                String AttributeType = "exon"
                String GTFAttributeType = "gene_id"
                String Stranded = "1"
        }
	prefix = basename(bam, ".bam")
        out = "${out_dir}"+"/"+"${prefix}"+".counts"
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
