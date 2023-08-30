version 1.0

import "../../wdl_tasks/trim_galore.wdl" as trimTasks
import "../../wdl_tasks/bowtie2.wdl" as bowtieTasks
import "../../wdl_tasks/bwa.wdl" as bwaTasks
import "../../wdl_tasks/filter_steps.wdl" as filterTasks
import "../../wdl_tasks/feature_count.wdl" as quantTasks
import "../../wdl_tasks/peak_calling.wdl" as pcTasks
import "../../wdl_tasks/make_bigWig.wdl" as makeBWWorkflow

workflow cut_and_run {
	input {
		File sampleList
		String fastq_dir
		String sample_out_dir
		String sampleName
		Boolean paired
		String Aligner 
		String? BWAIndex
		String PeakCaller
		String? PeakCallingControl
		Float? TopPeakFraction
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
	if ("${Aligner}" == "bwa") {
		call bwaTasks.BWA {
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
				sam=BWA.rawSam,
				sample_name="${sample_out_dir}"+"/"+"${sampleName}"+".raw"
		}
	}
	call filterTasks.remove_scaffolds {
		input:
			bam=sam_to_bam.bam,
			bam2=bowtie2.bam,
			chrom_no_scaff=ChromNoScaffold,
			sample_name=sampleName
	}
	call filterTasks.remove_blacklist {
		input:
			bam=remove_scaffolds.bam_noScaffold,
			blacklist=Blacklist,
			sample_name="${sample_out_dir}"+"/"+"${sampleName}"
	}
	call filterTasks.sort_bam {
		input:
			bam=remove_blacklist.bam_noBlacklist,
			sample_name="${sample_out_dir}"+"/"+"${sampleName}"
	}
	call filterTasks.remove_duplicates {
		input:
			bam=sort_bam.bam_sorted,
			sample_name="${sample_out_dir}"+"/"+"${sampleName}"
	}
	call filterTasks.size_filter_bam {
		input:
			bam=remove_duplicates.bam_noDuplicate,
			sample_name="${sample_out_dir}"+"/"+"${sampleName}",
			threshold_low=size_low,
			threshold_hi=size_high
	}
	if ("${PeakCaller}" == "macs2") {
		call pcTasks.macs2 {
			input:
				sampleName=sampleName,
				sample_out_dir=sample_out_dir,
				bam=size_filter_bam.hi,
				control_bam=PeakCallingControl
		}
		
	}
	if ("${PeakCaller}" == "seacr") {
		call pcTasks.bamToBedgraph {
			input:
				bam=size_filter_bam.hi,
				sampleName="${sample_out_dir}"+"/"+"${sampleName}"
		}
		call pcTasks.SEACR {
			input:
				sampleName=sampleName,
				sample_out_dir=sample_out_dir,
				bedgraph=bamToBedgraph.seacr_bg,
				control_bedgraph=PeakCallingControl,
				top_peak_fraction=TopPeakFraction
		}
	}
	call makeBWWorkflow.makeBigWig as BW_allReads {
		input:
			bam=remove_duplicates.bam_noDuplicate,
			chromosome_sizes=ChromosomeSizes,
			sampleName="${sample_out_dir}"+"/"+"${sampleName}.AllReads"
	}
	call makeBWWorkflow.makeBigWig as BW_hiThresh {
		input:
			bam=size_filter_bam.hi,
			chromosome_sizes=ChromosomeSizes,
			sampleName="${sample_out_dir}"+"/"+"${sampleName}.HighThreshold"
	}
	#call makeBWWorkflow.makeBigWig as BW_lowThresh {
        #        input:
        #                bam=size_filter_bam.low,
        #                chromosome_sizes=ChromosomeSizes,
        #                sampleName="${sample_out_dir}"+"/"+"${sampleName}.LowThreshold"
        #}	
	output {
		File bam = remove_duplicates.bam_noDuplicate
		File bw = BW_allReads.bw
		File? bam_low = size_filter_bam.low
		File? bam_hi = size_filter_bam.hi
		File? bw_hiThresh = BW_hiThresh.bw
		#File? bw_lowThresh = BW_lowThresh.bw
	}
}
