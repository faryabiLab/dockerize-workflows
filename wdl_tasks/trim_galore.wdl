version 1.0
# testing


# auto-zipped
task fastqc_trim {
	input {
		File? R1
		File? R2
		File? SE
		Boolean paired
		String sampleName
		String? Dockerhub_Pull = "faryabilab/trim_galore:0.10"
		Int quality = 15
		Int stringency = 5
		Float e = 0.1
		Int length = 20
		Boolean? illumina = True

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
			~{if illumina then "--illumina" else ""} \
			-e ${e} \
			--basename ${sampleName} \
			--length ${length} \
			${R1} ${R2}
		else
			trim_galore \
			-j ${cpu} \
			-q ${quality} \
			--fastqc \
			--phred33 \
			--gzip \
			--stringency ${stringency} \
			~{if illumina then "--illumina" else ""} \
			-e ${e} \
			--basename ${sampleName} \
			--length ${length} \
			${SE}
		fi
	}
	output {
		File? out_fqc  = "${sampleName}_trimmed.fq.gz"
    		File? out_fqc1 = "${sampleName}_val_1.fq.gz"
    		File? out_fqc2 = "${sampleName}_val_2.fq.gz"
    		File? out_html  = "${sampleName}_trimmed_fastqc.html"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
