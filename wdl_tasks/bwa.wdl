version 1.0

task BWA_paired {
	input {
		String sampleName
		String sample_out_dir

		String fastq1_trimmed
		String fastq2_trimmed
		String BWAIndex

		Int subsequence_seed = 32
		Int seed_max_edit_distance = 2
		Int read_trimming = 5
		
		Int MaximumInsertSize = 500
		Int MaximumReadOccurences = 100000
		Int MaxAlignmentsForXATag = 3
		Int MaxAlignmentsForXATag_DiscordanantPairs = 10
	}
	String readGroup = "@RG\tID:${sampleName}\tSM:${sampleName}"
	command {
		bwa aln \
		-q "${read_trimming}" \
		-l "${subsequence_seed}" \
		-k "${seed_max_edit_distance}" \
		"${BAWIndex}" \
		"${fastq1_trimmed}" > "${sample_out_dir}/"+sampleName+".1.sai"

		bwa aln \
                -q "${read_trimming}" \
                -l "${subsequence_seed}" \
                -k "${seed_max_edit_distance}" \
                "${BAWIndex}" \
                "${fastq2_trimmed}" > "${sample_out_dir}/"+sampleName+".2.sai"

		bwa sampe \
		-o "${MaximumReadOccurences}" \
		-a "${MaximumInsertSize}" \
		-n "${MaxAlignmentsForXATag}" \
		-N "${MaxAlignmentsForXATag_DiscordanantPairs}" \
		-r "${readGroup}" \
		"${sample_out_dir}/"+sampleName+".1.sai" \
		"${sample_out_dir}/"+sampleName+".2.sai" \
		> "${sample_out_dir}/"+"${sampleName}.raw.bam"

		rm "${sample_out_dir}/"+sampleName+".1.sai" "${sample_out_dir}/"+sampleName+".2.sai"
	}
	output {
		String rawBam = "${sample_out_dir}/"+"${sampleName}.raw.bam"
	}
}
