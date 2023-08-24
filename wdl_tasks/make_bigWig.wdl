version 1.0

####
# ALL arguments are REQUIRED
####

task read_count {
	input {
		String bam
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
		echo "scale=10; 1000000/${count}" | bc -l
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
		String bam
		String chromosome_sizes
		Float factor
		String sample_name
	}
	command {
		bamToBed -i "${bam}" -bed12 \
		| bed12ToBed6 -i stdin \
		| genomeCoverageBed -bg -i - -g ${chromosome_sizes} -scale ${factor} \
		| sort -k1,1 -k2,2n > "${sample_name}.bg"
	}
	output {
		File bedgraph = "${sample_name}.bg"
	}
	runtime {
		docker: 'faryabilab/bedtools:0.1.0'
	}
}

task bedgraph_to_bigwig {
	input {
                String bedgraph
		String chromosome_sizes
                String sample_name
        }
        command {
		bedGraphToBigWig ${bedgraph} ${chromosome_sizes} "${sample_name}.bw"
        }
        output {
                File bw = "${sample_name}.bw"
        }
        runtime {
                docker: 'faryabilab/bedtools:0.1.0'
        }
}
