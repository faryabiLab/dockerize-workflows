version 1.0

task read_count {
	input {
		File bam
	}	
	command {
		samtools view -c ${bam}
	}
	output {
		Int count = read_int(stdout())
	}
	runtime {
		docker: 'faryabilab/samtools:0.1.0'
	}
}

task calculate_factor {
	input {
		Int count
	}
	command {
		echo "scale=10; 1000000/${rc}" | bc -l
	}
	output {
		Float factor = read_float(stdout())
	}
	runtime {
		docker: 'faryabilab/bedtools:0.1.0'
	}
}

task bam_to_bedgraph {
	input {
		File bam
		File chromosome_sizes
		Float factor
		String out_dir
		String sample_name
	}
	out = "${out_dir}"+"/"+"${sample_name}"+".bg"
	command {
		bamToBed -i ${bam} -bed12 \
			| bed12ToBed6 -i stdin \
			| genomeCoverageBed -bg -i - -g ${chromosome_sizes} -scale ${factor} \
			| sort -k1,1 -k2,2n > ${out}
	}
	output {
		File bedgraph = ${out}
	}
	runtime {
		docker: 'faryabilab/bedtools:0.1.0'
	}
}

task bedgraph_to_bigwig {
	input {
                File bedgraph
		File chromosome_sizes
                String out_dir
                String sample_name
        }
        out = "${out_dir}"+"/"+"${sample_name}"+".bw"
        command {
		bedGraphToBigWig ${bedgraph} ${chromosome_sizes} ${out}
        }
        output {
                File bw = ${out}
        }
        runtime {
                docker: 'faryabilab/bedtools:0.1.0'
        }
}

}
