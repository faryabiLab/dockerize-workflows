version 1.0

import "../../wdl_tasks/trim_galore.wdl" as trimTasks
import "../../wdl_tasks/bwa.wdl" as bwaTasks
import "../../wdl_tasks/filter_steps.wdl" as filterTasks
import "../../wdl_tasks/peak_calling.wdl" as pcTasks
import "../../wdl_tasks/feature_count.wdl" as quantTasks
import "../../wdl_tasks/make_bigWig.wdl" as bwtasks

workflow chipseq {
	input {
		Boolean paired
                String fastq_dir
		String sample_out_dir
                String sampleName       
                String bwa_index
                String chromNoScaffold
		String blacklist
		String chromosome_sizes
		String? PeakcallingControl
	}
	call trimTasks.fastqc_trim {
		input:
			fastq_dir=fastq_dir,
			sample_out_dir=sample_out_dir,
			sampleName=sampleName,
			paired=paired
	}
	call bwaTasks.BWA {
		input:
			paired=paired,
			sampleName=sampleName,
			sample_out_dir=sample_out_dir,
			fastq1_trimmed=fastqc_trim.out_fqc1,
			fastq2_trimmed=fastqc_trim.out_fqc2,
			fastq_trimmed_single=fastqc_trim.out_fqc,
			BWAIndex=bwa_index
	}
	call filterTasks.sam_to_bam {
		input:
			sam=BWA.rawSam,
			sample_name="${sample_out_dir}/${sampleName}.raw"
	}
	call filterTasks.remove_scaffolds {
		input:
			bam=sam_to_bam.bam,
			chrom_no_scaff=chromNoScaffold,
			sample_name=sample_out_dir+"/"+sampleName
	}
	call filterTasks.remove_duplicates {
		input:
			bam=remove_scaffolds.bam_noScaffold,
			sample_name=sample_out_dir+"/"+sampleName,
	}
	call filterTasks.remove_blacklist {
		input:
			bam=remove_duplicates.bam_noDuplicate,
			blacklist=blacklist,
			sample_name=sample_out_dir+"/"+sampleName
	}
	call filterTasks.sort_bam {
		input:
			bam=remove_blacklist.bam_noBlacklist,
			sample_name=sample_out_dir+"/"+sampleName
	}
	call filterTasks.index_bam {
		input:
			bam=sort_bam.bam_sorted,
			sample_name=sample_out_dir+"/"+sampleName
	}
	call pcTasks.macs2 {
		sampleName=sampleName,
		sample_out_dir=sample_out_dir,
		bam=sort_bam.bam_sorted,
		control_bam=PeakcallingControl
	}
}
