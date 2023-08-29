version 1.0

import "../../wdl_tasks/trim_galore.wdl" as trimTasks
import "../../wdl_tasks/bowtie2.wdl" as bowtieTasks
import "../../wdl_tasks/bwa.wdl" as bwaTasks
import "../../wdl_tasks/filter_steps" as filterTasks
import "../../wdl_tasks/feature_count" as quantTasks
import "../../wdl_tasks/peak_calling.wdl" as pcTasks
import "../../wdl_tasks/make_bigWig.wdl" as bwTasks

workflow cut_and_run {
	input {
		File sampleList
		String fastq_dir
		String sample_out_dir
		String sampleName
		Boolean paired
		String Aligner 
		String? BWAIndex
		String? BowtieIndex
		String? PeakCaller
		String PeakCallingControl
		String ChromNoScaffold
		String Blacklist
		String ChromosomeSizes
		Int? size_low
		Int? size_high
	}
	call trimTasks.fastqc_trim {
		input:
			fastq_dir=fastq_dir,
			sample_out_dir=sample_out_dir,
			sampleName=sampleName,
			paired=paired
	}
	if (${Aligner} == "bowtie2") {
		call bowtieTasks.bowtie2 {
			input:
				sampleName=sampleName,
				sample_out_dir=sample_out_dir,
				fastq1_trimmed=fastqc_trim.out_fqc1,
				fastq2_trimmed=fastqc_trim.out_fqc2
		}
		call filterTasks.remove_scaffolds {
			input:
				bam=bowtie2.bam,
				chrom_no_scaff=ChromNoScaffold,
				sample_name=sampleName
		}
	}
	if (${Aligner} == "bwa") {
                call bwaTasks.bwa {
                        input:
				paired=paired,
				sampleName=sampleName,
				sample_out_dir=sample_out_dir,
				fastq1_trimmed=fastqc_trim.out_fqc1,
				fastq2_trimmed=fastqc_trim.out_fqc2,
				fastq_trimmed_single=fastqc_trim.out_fqc,
				BWAIndex=BWAIndex
                }
		call filterTasks.sam_to_bam {
			input:
				sam=bwa.rawSam,
				sample_name=${sample_out_dir}+"/"+${sampleName}.raw
		}
                call filterTasks.remove_scaffolds {
                        input:
                                bam=sam_to_bam.bam,
                                chrom_no_scaff=ChromNoScaffold,
                                sample_name=${sample_out_dir}+"/"+${sampleName}
                }
	
	}
	call filterTasks.remove_blacklist {
		input:
			bam=remove_scaffolds.bam_noScaffold,
			blacklist=Blacklist,
			sample_name=${sample_out_dir}+"/"+${sampleName}
	}
	call filterTasks.sort_bam {
		input:
			bam=remove_blacklist.bam_noBlacklist,
			sample_name=${sample_out_dir}+"/"+${sampleName}
	}
	call filterTasks.remove_duplicates {
		input:
			bam=sort_bam.bam_sorted,
			sample_name=${sample_out_dir}+"/"+${sampleName}
	}
	call filterTasks.size_filter_bam {
		input:
			bam=remove_duplicates.bam_noDuplicates,
			sample_name=${sample_out_dir}+"/"+${sampleName},
			threshold_low=size_low,
			threshold_hi=size_high
	}
	
	
}
