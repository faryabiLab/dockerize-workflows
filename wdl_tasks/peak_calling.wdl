version 1.0

task macs2_control {
	input {
		#### REQUIRED
		String sampleName
		String sample_out_dir
		String bam
		String control_bam		
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
		String BroadCutoff = 0.1			
		String BufferSize = 100000

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
	}
	command {
		macs2 callpeaks \
		-t ${bam} \
		${Pvalue} \
		${Qvalue} \
		-c ${control_bam} \
		-g ${GenomeSize} \
		--outdir ${sample_out_dir}
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
		${CallSummits} \
		${FixBimodal} \
		${NoModel} \
		${SaveFragPileup} \
		${TagSize} \
		${DownSample} \
		${NoLambda} \
		${BroadPeaks} \
		${BroadCutoff} \
		${CutoffAnalysis}
	}
	output {
		File narrowPeak = "${sampleOutDirWithPrefix}_peaks.narrowPeak"
		File summits = "${sampleOutDirWithPrefix}_summits.bed"
	}
	runtime {
		docker: 'faryabilab/macs2:0.1.0'
	}	
}

task macs2_noControl {
	input {
		#### REQUIRED
                String sampleName
                String sample_out_dir
                String bam
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
                String BroadCutoff = 0.1
                String BufferSize = 100000

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

        }
        command {
		macs2 callpeaks \
                -t ${bam} \
                ${Pvalue} \
                ${Qvalue} \
                -g ${GenomeSize} \
                --outdir ${sample_out_dir}
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
                ${CallSummits} \
                ${FixBimodal} \
                ${NoModel} \
                ${SaveFragPileup} \
                ${TagSize} \
                ${DownSample} \
                ${NoLambda} \
                ${BroadPeaks} \
                ${BroadCutoff} \
                ${CutoffAnalysis}
        }
        output {
		File narrowPeak = "${sampleOutDirWithPrefix}_peaks.narrowPeak"
                File summits = "${sampleOutDirWithPrefix}_summits.bed"
        }
        runtime {
                docker: 'faryabilab/macs2:0.1.0'
        }
}

task SEACR {
	input {
		String sampleName
                String sample_out_dir

                String bedgraph
		String control_bedgraph

		String Normalization = "norm"
		String RunMode = "stringent"

	}
	String sampleOutDirWithPrefix = "${sample_out_dir}/${sampleName}"
	command {
		bash "SEACR-1.3/SEACR_1.3.sh" \
		${bedgraph} \
		${control_bedgraph} \
		${Normalization} \
		${sampleOutDirWithPrefix}
	}
	output {
		File seacr_out = "${sampleOutDirWithPrefix}.${RunMode}.bed"
	}
	runtime {
		docker: 'faryabiLab/seacr:0.1.0'
	}
}
