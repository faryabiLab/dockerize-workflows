version 1.0

task macs2 {
	input {
		#### REQUIRED
		String sampleName
		String? bam
		String? control_bam		
		####
		String? Dockerhub_Pull = "faryabilab/macs2:0.1.0"		

		Int? pcBam
		String? FixBimodal
		String? NoModel
		Int? Shift
		Int? ExtensionSize
		String? SaveFragPileup
		String? TagSize
		String? Qvalue = 0.05
		String? Pvalue
		String? DownSample
		String? NoLambda
		String? BroadPeaks
		String? CutoffAnalysis
		String? CallSummits
		
		Int cpu = 12
		Int mem = 25
	}
	command {
		macs2 callpeak \
		-t ${bam} \
		~{if defined(Pvalue) then "-p "+ Pvalue else ""} \
		~{if defined(Qvalue) then "-q "+ Qvalue else ""} \
		~{if defined(control_bam) then "-c "+ control_bam else ""} \
		-n ${sampleName} \
		~{if defined(Shift) then "--shift "+ Shift else ""} \
		~{if defined(ExtensionSize) then "--extsize "+ ExtensionSize else ""} \
		--bw 300 \
		--seed 0 \
		--keep-dup 1 \
		${CallSummits} \
		${FixBimodal} \
		${NoModel} \
		${SaveFragPileup} \
		${TagSize} \
		${DownSample} \
		${NoLambda} \
		${BroadPeaks} \
		${CutoffAnalysis}
	}
	output {
		File narrowPeak = "${sampleName}_peaks.narrowPeak"
		File summits = "${sampleName}_summits.bed"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}	
}

task SEACR {
	input {
		String sampleName
                String? sample_out_dir
		String? Dockerhub_Pull = "faryabilab/seacr:0.1.0"

                String bedgraph
		String? control_bedgraph
		Float? top_peak_fraction

		String Normalization = "norm"
		String RunMode = "stringent"

		Int cpu = 12
		Int mem = 25

		String type

	}
	command {
		bash "/tmp/SEACR-1.3/SEACR_1.3.sh" \
		"${bedgraph}" \
		~{if defined(control_bedgraph) then control_bedgraph else top_peak_fraction} \
		"${Normalization}" \
		"${RunMode}" \
		"${sampleName}.${type}"
	}
	output {
		File seacr_out = "${sampleName}.${type}.${RunMode}.bed"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task bamToBedgraph {
	input {
		String? bam
		String sampleName
		String type
		String? Dockerhub_Pull = "faryabilab/bedtools:0.1.0"

		Int cpu = 12
		Int mem = 25
	}
	command {
		bedtools genomecov \
		-ibam "${bam}" \
		-bg \
		> "${sampleName}.seacr.bg.${type}"
	}
	output {
		File seacr_bg = "${sampleName}.seacr.bg.${type}"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
