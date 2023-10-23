version 1.0

task fastqc_trim {
	input {
		#### REQUIRED
		String fastq_dir
		String sample_out_dir
		String sampleName
		String? fastq_suffix
		Boolean paired
		####
		Int quality = 15
		Int stringency = 5
		Float e = 0.1
		Int length = 20

		Int cpu = 8
		Int mem = 8
	}
	command {
		if [[ "${paired}" == "true" ]]; then
			trim_galore \
			-j ${cpu} \
			-q ${quality} \
			--paired \
			--fastqc \
			--phred33 \
			--gzip \
			--stringency ${stringency} \
			-e ${e} \
			--basename ${sampleName} \
			--length ${length} \
			-o ${sample_out_dir} \
			"${fastq_dir}/${sampleName}_R1${fastq_suffix}.fastq.gz" \
			"${fastq_dir}/${sampleName}_R2${fastq_suffix}.fastq.gz"
		else
			trim_galore \
			-j ${cpu} \
			-q ${quality} \
			--fastqc \
			--phred33 \
			--gzip \
			--stringency ${stringency} \
			-e ${e} \
			--basename ${sampleName} \
			--length ${length} \
			-o ${sample_out_dir} \
			"${fastq_dir}/${sampleName}.fastq.gz"
		fi
	}
	output {
		String? out_fqc = "${sample_out_dir}/${sampleName}"+"_trimmed.fq.gz"
		String? out_fqc1 = "${sample_out_dir}/${sampleName}"+"_val_1.fq.gz"
		String? out_fqc2 = "${sample_out_dir}/${sampleName}"+"_val_2.fq.gz"
	}
	runtime {
		docker: "faryabilab/trim_galore:0.10"
		cpu: 8
		mem: 16
	}
}
