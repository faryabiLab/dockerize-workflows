version 1.0

##########
# Align RNAseq fastq files with STAR
##########
task STAR {
	input {
		#### REQUIRED
		File? fastq1_trimmed
		File? fastq2_trimmed			
		File? fastq_trimmed_single
		File star_index
		String sample_name
		Boolean paired
		####
		String? Dockerhub_Pull = "faryabilab/star:0.10"
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
		Int STAR_cpu=10
		Int STAR_mem=25

	}
	command {
		
		mkdir star_index
        	tar -xzf ~{star_index_tar} -C star_index
		
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
	output {File bam = "${sample_name}Aligned.out.bam"}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${STAR_cpu}"
		mem: "${STAR_mem}"
	}
}

##########
# Zip STAR index files for easy localization
##########
task prep_star_index
{
	input
	{
		String star_index_dir
	}
	command
	{
		tar -czf ${star_index_dir}
	}
	output 
	{
		File index_tar = "star_index.tar.gz"
	}
	runtime 
	{
		docker: "ubuntu:20.04"
    		cpu: 1
    		memory: "4G"
	}
}


