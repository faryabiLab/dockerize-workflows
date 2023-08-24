version 1.0 

task STAR {
	input {
		#### REQUIRED
		String fastq1_trimmed
		String? fastq2_trimmed			
		String star_index
		String sample_name
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
	}
	command {
		STAR \
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
	}
	output {File bam = "${sample_name}Aligned.out.bam"}
	runtime {
		docker: 'faryabilab/star:0.10'
	}
}
