version 1.0

task bowtie2 {
	input {
		String sampleName
		String sample_out_dir
		String fastq1_trimmed
		String fastq2_trimmed
		String BowtieIndex
		
		Int Seed = 0

		Int Trim5 = 0
		Int Trim3 = 0
		Int MaxMismatch = 0
		Int SeedSubstringLength = 22
		Int ExtraRefChars = 15
		Int DisallowGaps = 4
		Int MaxFragmentLength = 500
		Int MinFragmentLength = 0
		
		String? IgnoreQualities
		String? NoForwardAlign
		String? NoReverseAlign
		String? EndtoEnd = "--end-to-end"
		String? LocalAlignment
		String? NoMixed
		String? NoDiscordant
		String? DoveTail
		String? NoContain
		String? NoOverlap

		#Output args
		String MetricsFile = "bowtie2MetricsFile.txt"
	}
	sampleOutDirWithPrefix = ${sample_out_dir}+"/"+${sampleName}

	String IgnoreQualities = ${default="" IgnoreQualities}
	String NoForwardAlign = ${default="" NoForwardAlign}
	String NoReverseAlign = ${default="" NoReverseAlign}
	String EndtoEnd = ${default="--end-to-end" EndtoEnd}
	String LocalAlignment = ${default="" LocalAlignment}
	String NoMixed = ${default="" NoMixed}
	String NoDiscordant ${default="" NoDiscordant}
	String DoveTail = ${default="" DoveTail}
	String NoContain = ${default="" NoContain}
	String NoOverlap = ${default = "" NoOverlap}
	String MetricsFile = "${sampleoutDirWithPrefix}"+"."+"${MetricsFile}"	
	command {
		bowtie2 \
		-x ${BowtieIndex} \
		-1 ${fastq1_trimmed} -2 ${fastq2_trimmed} \
		-S "${sampleOutDirWithPrefix}.out.tmp.sam" \
		--trim5 ${Trim5} --trim3 ${Trim3} \
		-N ${MaxMismatch} \
		-L ${SeedSubstringLength} \
		--dpad ${ExtraRefChars} \
		--gbar ${DisallowGaps} \
		--met-file ${MetricsFile} \
		${EndtoEnd} \
		${IgnoreQualities} \
		${NoForwardAlign} \
		${NoReverseAlign} \
		${LocalAlignment} \
		${NoMixed} \
		${NoDiscordant} \
		${DoveTail} \
		${NoContain} \
		${NoOverlap} 
		samtols view -b "${sampleOutDirWithPrefix}.out.tmp.sam" > "${sampleOutDirWithPrefix}.raw.bam"
	
		rm "${sampleOutDirWithPrefix}.out.tmp.sam"	
	}
	output {
		File bam = "${sampleOutDirWithPrefix}.raw.bam"
	}
	runtime {
		docker: 'faryabilab/bowtie2:0.1.1'
	}
}
