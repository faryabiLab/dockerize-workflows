version 1.0 

task STAR {
	input {
		File fastq1_trimmed
		File? fastq2_trimmed			
		String star_index
		String out_dir
		String sample_name
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
	String out="${out_dir}"+"/"+"${sample_name}"
	command {
		STAR \
		--genomeDir ${star_index} \
		--outFileNamePrefix ${out} \
		--outFilterType ${outFilterType} \
		--readFilesIn ${fastq1_trimmed} ${fastq2_trimmed} \
		--readFilesCommand ${readFilesCommand} \
		--outFilterIntronMotifs ${outFilterIntronMotifs} \
		--alignIntronMax ${alignIntronMax} \
		--outSamstrandField ${outSamstrandField} \
		--outSAMunmapped ${outSAMunmapped} \
		--chimSegmentMin ${cimSegmentMin} \
		--chimJunctionOverhangMin ${chimJunctionOverhangMin} \
		--outSAMtype ${outSAMtype}
	}
	output {
		File bam = "${out}.out.bam"
	}
	runtime {
		docker: 'faryabilab/star:0.10'
	}
}
