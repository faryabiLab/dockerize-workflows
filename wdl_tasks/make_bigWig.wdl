version 1.2

####
# ALL arguments are REQUIRED
####

task read_count {
	input {
		File? bam
		String? Dockerhub_Pull = "faryabilab/samtools:0.1.0"	
		Int cpu = 1
		Int mem = 5
	}	
	command {
		samtools view -c ${bam}
	}
	output {
		Int count = read_int(stdout())
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task calculate_factor {
	input {
		Int count
		String? Dockerhub_Pull = "faryabilab/bedtools:0.1.0"

		Int cpu = 1
		Int mem = 5
	}
	command {
		echo "scale=10; 1000000/${count}" | bc -l
	}
	output {
		Float factor = read_float(stdout())
	}
	runtime {
		docker: "${Dockerhub_Pull}"
                cpu: "${cpu}"
                mem: "${mem}"
	}
}

task bam_to_bedgraph {
	input {
		File? bam
		String chromosome_sizes
		Float factor
		String sample_name
		String? Dockerhub_Pull = "faryabilab/bedtools:0.1.0"
		
		Int cpu = 8
		Int mem = 16
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
		docker: "${Dockerhub_Pull}"
                cpu: "${cpu}"
                mem: "${mem}"
	}
}

task bedgraph_to_bigwig {
	input {
                File bedgraph
		String chromosome_sizes
                String sample_name
		String? Dockerhub_Pull = "faryabilab/bedtools:0.1.0"
		
		Int cpu = 12
		Int mem = 16
        }
        command {
		bedGraphToBigWig ${bedgraph} ${chromosome_sizes} "${sample_name}.bw"
		rm ${bedgraph}
        }
        output {
                File bw = "${sample_name}.bw"
        }
        runtime {
                docker: "${Dockerhub_Pull}"
                cpu: "${cpu}"
                mem: "${mem}"
        }
}

workflow makeBigWig {
	input {
		File? bam
		String chromosome_sizes
		String sampleName
	}
	call read_count {
		input:
			bam=bam
	}
	call calculate_factor {
		input:
			count=read_count.count
	}
	call bam_to_bedgraph {
		input:
			bam=bam,
			chromosome_sizes=chromosome_sizes,
			factor=calculate_factor.factor,
			sample_name=sampleName
	}
	call bedgraph_to_bigwig {
		input:
			bedgraph=bam_to_bedgraph.bedgraph,
			chromosome_sizes=chromosome_sizes,
			sample_name=sampleName
	}
	output {
		File bw = bedgraph_to_bigwig.bw
      
	}
}
