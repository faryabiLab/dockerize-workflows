version 1.0
task BWA {
	input {
		#### REQUIRED
		Boolean paired
		String sampleName
		String? fastq1_trimmed
		String? fastq2_trimmed
		String? fastq_trimmed_single
		String BWAIndex
		####
		
		String? Dockerhub_Pull="faryabilab/bwa:0.1.0"

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
		Int? DiscardMemOccurences = 10000
		Int? MatchingScore = 1
		Int? MismatchPenalty = 4		
		Int? GapOpenPenalty = 6
		Int? GapExtensionPenalty = 1
		Int? ClippingPenality = 5
		Int? ScoreCutoff = 30
		String? HardClipping
		String? OutputUnpairedReads

		Int cpu = 16
		Int mem = 32	
	}
	command {
		if [[ "${paired}" == "true" ]]; then 
			bwa mem \
			-t ${cpu} \
			-k ${MinSeedLength} \
			-w ${Bandwidth} \
			-d ${ZDropoff} \
			-r ${TriggerReSeed} \
			-c ${DiscardMemOccurences} \
			-A ${MatchingScore} \
			-B ${MismatchPenalty} \
			-O ${GapOpenPenalty} \
			-E ${GapExtensionPenalty} \
			-L ${ClippingPenality} \
			-R "@RG\tID:${sampleName}\tSM:${sampleName}" \
			-T ${ScoreCutoff} \
			~{if defined(HardClipping) then "-H" else ""} \
			~{if defined(OutputUnpairedReads) then "-a" else ""} \
			"${BWAIndex}" \
			"${fastq1_trimmed}" "${fastq2_trimmed}" \
			> "${sampleName}.raw.sam"	


		else
			bwa mem \
			-t ${cpu} \
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
			> "${sampleName}.raw.sam"

		fi
			
	}
	output {
		File rawSam = "${sampleName}.raw.sam"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
