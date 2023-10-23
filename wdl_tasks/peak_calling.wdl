version 1.0

task macs2 {
	input {
		#### REQUIRED
		String sampleName
		String sample_out_dir
		String? bam
		String? control_bam		
		####		

		String GenomeSize = 'hs'
		Int Shift = 0		
		Int ExtensionSize = 20
		Int BandWidth = 200
		Int MinFragSize = 20
		String MFold = "5 50"		
		Int Seed = 0 
                Int SmallLocal = 1000
                Int LargeLocal = 10000		
		String ScaleTo = "small"
		Float BroadCutoff = 0.1			
		String BufferSize = 100000

		Int? pcBam
		String? FixBimodal
		String? NoModel
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
		-t  ${bam} \
		~{if defined(Pvalue) then "-p "+ Pvalue else ""} \
		~{if defined(Qvalue) then "-q "+ Qvalue else ""} \
		~{if defined(control_bam) then "-c "+ control_bam else ""} \
		-g ${GenomeSize} \
		--outdir ${sample_out_dir} \
		-n ${sampleName} \
		--shift ${Shift} \
		--extsize ${ExtensionSize} \
		--bw ${BandWidth} \
		--d-min ${MinFragSize} \
		--mfold ${MFold} \
		--seed ${Seed} \
		--slocal ${SmallLocal} \
		--llocal ${LargeLocal} \
		--scale-to ${ScaleTo} \
		--buffer-size ${BufferSize} \
		--broad-cutoff ${BroadCutoff} \
		--outdir ${sample_out_dir} \
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
		String narrowPeak = "${sample_out_dir}/${sampleName}_peaks.narrowPeak"
		String summits = "${sample_out_dir}/${sampleName}_summits.bed"
	}
	runtime {
		docker: 'faryabilab/macs2:0.1.0'
		cpu: "${cpu}"
		mem: "${mem}"
	}	
}

task SEACR {
	input {
		String sampleName
                String sample_out_dir

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
		"${bedgraph}.${type}" \
		~{if defined(control_bedgraph) then control_bedgraph else top_peak_fraction} \
		"${Normalization}" \
		"${RunMode}" \
		"${sample_out_dir}/${sampleName}.${type}"
	}
	output {
		String seacr_out = "${sample_out_dir}/${sampleName}.${type}.${RunMode}.bed"
	}
	runtime {
		docker: 'faryabilab/seacr:0.1.0'
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task bamToBedgraph {
	input {
		String? bam
		String sampleName
		String type

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
		String seacr_bg = "${sampleName}.seacr.bg"
	}
	runtime {
		docker: "faryabilab/bedtools:0.1.0"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
