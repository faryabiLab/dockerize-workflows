version 1.0
# auto-zipped
task fastqc_trim {
	input {
		File? R1
		File? R2
		File? SE
		Boolean paired
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
			"$R1" "$R2"
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
			"$SE"
		fi
	}
	output {
		Array[File] reports = glob("*.html")
		Array[File] trimmed = glob("*.fq.gz")
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
