version 1.0
## TESTT

import "wdl_tasks/trim_galore.wdl" as trimTasks
import "wdl_tasks/star.wdl" as starTasks
import "wdl_tasks/filter_steps.wdl" as filterTasks
import "wdl_tasks/feature_count.wdl" as quantifyTasks
import "wdl_tasks/make_bigWig.wdl" as bwTasks

workflow RNAseq 
{
	input 
	{
		String sampleList
		Boolean paired       
                File star_index
                File chromNoScaffold
		File GeneAnnotationFile
		File chromosome_sizes
		File? blacklist
		String? fastq_suffix
	}
	
	Array[Array[String]] samples = read_tsv(sampleList)
	Boolean blDefined = defined(blacklist)

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
	
		call starTasks.STAR 
		{
			input:
				fastq_trimmed_single =  trimmed_single,
				fastq1_trimmed = 	trimmed_r1,
				fastq2_trimmed =	trimmed_r2,
				star_index =		star_index,
				sample_name =		sample_id,
				paired =		paired
		}

		call filterTasks.remove_scaffolds 
		{
			input:
				bam=STAR.bam,
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
		if (paired) 
		{
			call quantifyTasks.count_reads_paired 
			{
				input:
					bam=sort_bam.bam_sorted,
					GeneAnnotationFile=GeneAnnotationFile,
					sample_name=sample_id
			}
		}
		if (!paired) 
		{
			call quantifyTasks.count_reads_single 
			{
				input:
					bam=sort_bam.bam_sorted,
					GeneAnnotationFile=GeneAnnotationFile,
					sample_name=sample_id
			}
		}
		call bwTasks.read_count 
		{
			input:
				bam=sort_bam.bam_sorted
		}
		call bwTasks.calculate_factor 
		{
			input:
				count=read_count.count
		}
		call bwTasks.bam_to_bedgraph 
		{
			input:
				bam=sort_bam.bam_sorted,
				chromosome_sizes=chromosome_sizes,
				factor=calculate_factor.factor,
				sample_name=sample_id
		}
		call bwTasks.bedgraph_to_bigwig 
		{
			input:
				bedgraph=bam_to_bedgraph.bedgraph,
				chromosome_sizes=chromosome_sizes,
				sample_name=sample_id
		}	
        }	
	output 
	{
		Array[File] finalBam = sort_bam.bam_sorted
		Array[File] finalBamIndex = index_bam.bam_index
		Array[File?] countsSingle = count_reads_single.gene_counts
		Array[File?] countsPaired = count_reads_paired.gene_counts
		Array[File] bw = bedgraph_to_bigwig.bw
	}
}
