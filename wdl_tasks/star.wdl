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
		
		String outFilterType
		String readFilesCommand
		String outSamAttributes
		String outFilterIntronMotifs
		Int alignIntronMax
		String outSamstrandField
		String outSAMunmapped
		Int chimSegmentMin
		Int chimJunctionOverhangMin
		String outSAMtype

		Int STAR_cpu
		Int STAR_mem
	}
	command {
		if [[ "~{paired}" == "true" ]]; then
			STAR \
			--runThreadN ${STAR_cpu} \
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
			--runThreadN ${STAR_cpu} \
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
	output {String bam = "${sample_name}Aligned.out.bam"}
	runtime {
		docker: 'faryabilab/star:0.10'
		cpu: 10
		mem: 24
	}
}
