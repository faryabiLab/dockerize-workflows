version 1.0
task BWA {
	input {
		#### REQUIRED
		Boolean paired
		String sampleName
		String sample_out_dir
		String? fastq1_trimmed
		String? fastq2_trimmed
		String? fastq_trimmed_single
		String BWAIndex
		####

		Int? subsequence_seed = 32
		Int? seed_max_edit_distance = 2
		Int? read_trimming = 5
		
		Int? MaximumInsertSize = 500
		Int? MaximumReadOccurences = 100000
		Int? MaxAlignmentsForXATag = 3
		Int? MaxAlignmentsForXATag_DiscordanantPairs = 10

		Int? MinSeedLength = 19
		Int? Bandwidth = 100
		Int? ZDropoff = 100
		Float? TriggerReSeed = 1.5
		Int? MatchingScore = 1
		Int? MismatchPenalty = 4		
		Int? GapOpenPenalty = 6
		Int? GapExtensionPenalty = 1
		Int? ClippingPenality = 5
		Int? ScoreCutoff = 30
		String? HardClipping	
	}
	command {
		if [[ "${paired}" == "true" ]]; then 
			bwa mem \
			-k ${MinSeedLength} \
			-w ${Bandwidth} \
			-d ${ZDropoff} \
			-r ${TriggerReSeed} \
			-A ${MatchingScore} \
			-B ${MismatchPenalty} \
			-O ${GapOpenPenalty} \
			-E ${GapExtensionPenalty} \
			-L ${ClippingPenality} \
			-R "@RG\tID:${sampleName}\tSM:${sampleName}" \
			-T ${ScoreCutoff} \
			"${BWAIndex}" \
			"${fastq1_trimmed}" "${fastq2_trimmed}" \
			> "${sample_out_dir}/${sampleName}.raw.sam"	


		else
			bwa mem \
			-k ${MinSeedLength} \
			-w ${Bandwidth} \
			-d ${ZDropoff} \
			-r ${TriggerReSeed} \
			-A ${MatchingScore} \
			-B ${MismatchPenalty} \
			-O ${GapOpenPenalty} \
			-E ${GapExtensionPenalty} \
			-L ${ClippingPenality} \
			-R "@RG\tID:${sampleName}\tSM:${sampleName}" \
			-T ${ScoreCutoff} \
			"${BWAIndex}" \
			"${fastq_trimmed_single}" \
			> "${sample_out_dir}/${sampleName}.raw.sam"

		fi
			
	}
	output {
		String rawSam = "${sample_out_dir}/${sampleName}.raw.sam"
	}
	runtime {
		docker: 'faryabilab/bwa:0.1.0'
		cpu: 10
		memory: 10
	}
}
