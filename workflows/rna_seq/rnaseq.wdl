version 1.0

import "wdl_tasks/trim_galore.wdl" as trimTasks
import "wdl_tasks/star.wdl" as starTasks
import "wdl_tasks/filter_steps.wdl" as filterTasks
import "wdl_tasks/feature_count.wdl" as quantifyTasks
import "wdl_tasks/make_bigWig.wdl" as bwTasks

#####
# The model for output file naming is as follows:
# Each sample name will be a directory, and wthin each directory will be that sample's complete collection of files
# from the entire run of the pipeline
#####

workflow RNAseq 
{
	input 
	{
		String sampleList
		String fastq_dir
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
		String sampleName=sample[0]
		
		# Build paths as Strings (no optionals, no None/null here)
   		File? R1a_path = fastq_dir + "/" + sampleName + "_R1" + select_first([fastq_suffix, ""]) + ".fastq.gz"
    		File? R2a_path = fastq_dir + "/" + sampleName + "_R2" + select_first([fastq_suffix, ""]) + ".fastq.gz"
    		File? R1b_path = fastq_dir + "/" + sampleName + "_1"  + select_first([fastq_suffix, ""]) + ".fastq.gz"
    		File? R2b_path = fastq_dir + "/" + sampleName + "_2"  + select_first([fastq_suffix, ""]) + ".fastq.gz"
    		File? SE_path  = fastq_dir + "/" + sampleName + ".fastq.gz"

    		# Choose naming convention with select_first; still no None/null
    		File? R1_choice = select_first([R1a_path, R1b_path])
    		File? R2_choice = select_first([R2a_path, R2b_path])
    		File? SE_choice = SE_path
		
		if (paired) 
		{
     			call trimTasks.fastqc_trim as trim_pe
			{
       				input:
          				R1 = R1_choice,
          				R2 = R2_choice,
          				paired = paired,
          				sampleName = sampleName
      			}
    		}

    		if (!paired) 
		{
      			call trimTasks.fastqc_trim as trim_se
			{
        			input:
          				SE = SE_choice,
          				paired = paired,
          				sampleName = sampleName
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
				sample_name =		sampleName,
				paired =		paired
		}

		call filterTasks.remove_scaffolds 
		{
			input:
				bam=STAR.bam,
				chrom_no_scaff=chromNoScaffold,
				sample_name=sampleName
		}
		call filterTasks.remove_duplicates 
		{
			input:
				bam=remove_scaffolds.bam_noScaffold,
				sample_name=sampleName
		}
		# Conditionally call remove_blacklist if blDefined = true
		if (blDefined)
		{
			call filterTasks.remove_blacklist
                	{
                        	input:
                                	bam=remove_duplicates.bam_noDuplicate,
                                	blacklist=blacklist,
                                	sample_name=sampleName
                	}

		}
		call filterTasks.sort_bam 
		{
			input:
				bam=select_first([remove_blacklist.bam_noBlacklist, remove_scaffolds.bam_noScaffold]),
				sample_name=sampleName
		}
		call filterTasks.index_bam 
		{
			input:
				bam=sort_bam.bam_sorted,
				sample_name=sampleName
		}
		if (paired) 
		{
			call quantifyTasks.count_reads_paired 
			{
				input:
					bam=sort_bam.bam_sorted,
					GeneAnnotationFile=GeneAnnotationFile,
					sample_name=sampleName
			}
		}
		if (!paired) 
		{
			call quantifyTasks.count_reads_single 
			{
				input:
					bam=sort_bam.bam_sorted,
					GeneAnnotationFile=GeneAnnotationFile,
					sample_name=sampleName
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
				sample_name=sampleName
		}
		call bwTasks.bedgraph_to_bigwig 
		{
			input:
				bedgraph=bam_to_bedgraph.bedgraph,
				chromosome_sizes=chromosome_sizes,
				sample_name=sampleName
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
