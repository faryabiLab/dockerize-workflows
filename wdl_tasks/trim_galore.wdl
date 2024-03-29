version 1.0

task fastqc_trim {
	input {
		#### REQUIRED
		String fastq_dir
		String sampleName
		String? fastq_suffix
		Boolean paired
		####
		String? Dockerhub_Pull = "faryabilab/trim_galore:0.10"
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
			"${fastq_dir}/${sampleName}.fastq.gz"
		fi
	}
	output {
		File? out_html = "${sampleName}"+"_trimmed_fastqc.html"
		File? out_fqc = "${sampleName}"+"_trimmed.fq.gz"
		File? out_fqc1 = "${sampleName}"+"_val_1.fq.gz"
		File? out_fqc2 = "${sampleName}"+"_val_2.fq.gz"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
