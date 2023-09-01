version 1.0 

task STAR {
	input {
		#### REQUIRED
		String? fastq1_trimmed
		String? fastq2_trimmed			
		String? fastq_trimmed_single
		String star_index
		String sample_name
		Boolean paired
		####
		
		String outFilterType = "BySJout"
		String readFilesCommand = "zcat"
		String outSamAttributes = "Standard"
		String outFilterIntronMotifs = "RemoveNoncanonicalUnannotated"
		Int alignIntronMax = 100000
		String outSamstrandField = "intronMotif"	
		String outSAMunmapped = "Within"
		Int chimSegmentMin = 25
		Int chimJunctionOverhangMin = 25
		String outSAMtype = "BAM Unsorted"

		Int cpu = 16
		Int mem = 50
	}
	command {
		if [[ "~{paired}" == "true" ]]; then
			STAR \
			--runThreadN ${cpu} \
			--genomeDir ${star_index} \
			--readFilesIn ${fastq1_trimmed} ${fastq2_trimmed} \
			--outFileNamePrefix ${sample_name} \
			--readFilesCommand ${readFilesCommand} \
			--outSAMattributes ${outSamAttributes} \
			--outFilterIntronMotifs ${outFilterIntronMotifs} \
			--outFilterType ${outFilterType} \
			--alignIntronMax ${alignIntronMax} \
			--outSAMstrandField ${outSamstrandField} \
			--outSAMunmapped ${outSAMunmapped} \
			--chimSegmentMin ${chimSegmentMin} \
			--outSAMtype ${outSAMtype}
		else
			STAR \
			--runThreadN ${cpu} \
			--genomeDir ${star_index} \
			--readFilesIn ${fastq_trimmed_single} \
			--outFileNamePrefix ${sample_name} \
			--readFilesCommand ${readFilesCommand} \
			--outSAMattributes ${outSamAttributes} \
			--outFilterIntronMotifs ${outFilterIntronMotifs} \
			--outFilterType ${outFilterType} \
			--alignIntronMax ${alignIntronMax} \
			--outSAMstrandField ${outSamstrandField} \
			--outSAMunmapped ${outSAMunmapped} \
			--chimSegmentMin ${chimSegmentMin} \
			--outSAMtype ${outSAMtype}
		fi
	}
	output {File bam = "${sample_name}Aligned.out.bam"}
	runtime {
		docker: 'faryabilab/star:0.10'
		cpu: "${cpu}"
		#memory: "${mem}"
	}
}
