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

		Int? cpu = 8
	}
	command {
		if [[ "~{paired} == "true" ]]; then 
			bwa aln \
			-q "${read_trimming}" \
			-l "${subsequence_seed}" \
			-k "${seed_max_edit_distance}" \
			"${BWAIndex}" \
			"${fastq1_trimmed}" > "${sample_out_dir}/"+sampleName+".1.sai"

			bwa aln \
			-q "${read_trimming}" \
			-l "${subsequence_seed}" \
			-k "${seed_max_edit_distance}" \
			"${BWAIndex}" \
			"${fastq2_trimmed}" > "${sample_out_dir}/"+sampleName+".2.sai"

			bwa sampe \
			-o "${MaximumReadOccurences}" \
			-a "${MaximumInsertSize}" \
			-n "${MaxAlignmentsForXATag}" \
			-N "${MaxAlignmentsForXATag_DiscordanantPairs}" \
			-r "@RG\tID:${sampleName}\tSM:${sampleName}" \
			"${sample_out_dir}/"+sampleName+".1.sai" \
			"${sample_out_dir}/"+sampleName+".2.sai" \
			"${fastq1_trimmed}" "${fastq2_trimmed}" \
			> "${sample_out_dir}/"+"${sampleName}.raw.sam"

			rm "${sample_out_dir}/"+sampleName+".1.sai" "${sample_out_dir}/"+sampleName+".2.sai"
		else
			bwa mem \
			-t ${cpu} \
			-k ${MinSeedLength}
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
			
	}
	output {
		String rawBam = "${sample_out_dir}/"+"${sampleName}.raw.sam"
	}
	runtime {
		docker: 'faryabilab/bwa:0.1.0'
	}
}
