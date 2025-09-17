version 1.0

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
		# Decide naming convention of fastqs

		if [[ "${paired}" == "true" ]]; then
			if [[ -f "${fastq_dir}/${sampleName}_R1${fastq_suffix}.fastq.gz" ]]; then
                                R1="${fastq_dir}/${sampleName}_R1${fastq_suffix}.fastq.gz"
                                R2="${fastq_dir}/${sampleName}_R2${fastq_suffix}.fastq.gz"
                        elif [[ -f "${fastq_dir}/${sampleName}_1${fastq_suffix}.fastq.gz" ]]; then
                                R1="${fastq_dir}/${sampleName}_1${fastq_suffix}.fastq.gz"
                                R2="${fastq_dir}/${sampleName}_2${fastq_suffix}.fastq.gz"
                        else
                                echo "ERROR: Cannot find FASTQ files for sample ${sampleName} with either _R1/_R2 or _1/_2 pattern" >&2
                                exit 1
			fi

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
