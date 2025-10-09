version 1.0

import "wdl_tasks/trim_galore.wdl" as trimTasks
import "wdl_tasks/bwa.wdl" as bwaTasks
import "wdl_tasks/filter_steps.wdl" as filterTasks
import "wdl_tasks/peak_calling.wdl" as pcTasks
import "wdl_tasks/feature_count.wdl" as quantTasks
import "wdl_tasks/make_bigWig.wdl" as bwTasks

workflow ChIPseq 
{
	input 
	{
		File sampleList
		Boolean paired
		File bwa_index
		File chromNoScaffold
		File chromosome_sizes
		File? blacklist
		File? PeakcallingControl
	}

	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) 
	{
		String sample_id = sample[0]
		File R1 = sample[1]
		File? R2 = if length(sample) > 2 && sample[2] != "" then sample[2] else ""
	
		if (paired) 
		{
     			call trimTasks.fastqc_trim as trim_pe
			{
       				input:
          				R1 = R1,
          				R2 = R2,
          				paired = paired,
          				sampleName = sample_id
      			}
    		}

    		if (!paired) 
		{
      			call trimTasks.fastqc_trim as trim_se
			{
        			input:
          				SE = R1,
          				paired = paired,
          				sampleName = sample_id
      			}
    		}

		File? trimmed_single = trim_se.out_fqc
		File? trimmed_r1     = trim_pe.out_fqc1
		File? trimmed_r2     = trim_pe.out_fqc2
	
		call bwaTasks.BWA 
		{
			input:
				paired=paired,
				sampleName=sample_id,
				fastq1_trimmed=trimmed_r1,
				fastq2_trimmed=trimmed_r2,
				fastq_trimmed_single=trimmed_single,
				BWAIndex=bwa_index,
		}
		call filterTasks.sam_to_bam 
		{
			input:
				sam=BWA.rawSam,
				sample_name="${sample_id}.raw"
		}
		call filterTasks.remove_unmapped
		{
			input:
				bam=sam_to_bam.bam,
				sample_name=sample_id
		}
		call filterTasks.remove_lowQuality
		{
			input:
				bam=remove_unmapped.bam_noUnmapped,
				sample_name=sample_id
		}
		call filterTasks.sort_bam_name
		{
			input:
				bam=remove_lowQuality.bam_noLowQuality,
				sample_name=sample_id
		}
		call filterTasks.fix_mate
		{
			input:
				bam=sort_bam_name.bam_sortedByname,
				sample_name=sample_id
				
		}
		call filterTasks.remove_scaffolds 
		{
			input:
				bam=fix_mate.bam_fixMate,
				chrom_no_scaff=chromNoScaffold,
				sample_name=sample_id
		}
		call filterTasks.mark_duplicates 
		{
			input:
				bam=remove_scaffolds.bam_noScaffold,
				sample_name=sample_id
		}
		call filterTasks.remove_duplicates_unmapped
		{
			input:
				bam=mark_duplicates.bam_dupMarked,
				sample_name=sample_id
		}
		
		Boolean blDefined = defined(blacklist)
		# Conditionally call remove_blacklist if blDefined = true
		if (blDefined)
		{
			call filterTasks.remove_blacklist
                	{
                        	input:
                                	bam=remove_duplicates_unmapped.bam_noDuplicate,
                                	blacklist=blacklist,
                                	sample_name=sample_id
                	}

		}
		call filterTasks.sort_bam 
		{
			input:
				bam=select_first([remove_blacklist.bam_noBlacklist, remove_duplicates_unmapped.bam_noDuplicate]),
				sample_name=sample_id
		}
		call filterTasks.index_bam 
		{
			input:
				bam=sort_bam.bam_sorted,
				sample_name=sample_id
		}
		call pcTasks.macs2 {
			input:
				sampleName=sample_id,
				bam=sort_bam.bam_sorted,
				control_bam=PeakcallingControl	
		}
		call bwTasks.read_count {
			input:
				bam=sort_bam.bam_sorted
		}
		call bwTasks.calculate_factor {
			input:
				count=read_count.count
		}
		call bwTasks.bam_to_bedgraph {
			input:
				bam=sort_bam.bam_sorted,
				chromosome_sizes=chromosome_sizes,
				factor=calculate_factor.factor,
				sample_name=sample_id
		}
		call bwTasks.bedgraph_to_bigwig {
			input:
				bedgraph=bam_to_bedgraph.bedgraph,
				chromosome_sizes=chromosome_sizes,
				sample_name=sample_id
		}
	}
	output {
		Array[File] finalBam = sort_bam.bam_sorted
		Array[File] finalBamIndex = index_bam.bam_index
		Array[File] peaks = macs2.narrowPeak
		Array[File] bw = bedgraph_to_bigwig.bw
	}
}
