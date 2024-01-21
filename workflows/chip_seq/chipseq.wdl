version 1.0

import "wdl_tasks/trim_galore.wdl" as trimTasks
import "wdl_tasks/bwa.wdl" as bwaTasks
import "wdl_tasks/filter_steps.wdl" as filterTasks
import "wdl_tasks/peak_calling.wdl" as pcTasks
import "wdl_tasks/feature_count.wdl" as quantTasks
import "wdl_tasks/make_bigWig.wdl" as bwTasks

#####
# The model for output file naming is as follows:
# Each sample name will be a directory, and wthin each directory will be that sample's complete collection of files
# from the entire run of the pipeline
#####

workflow ChIPseq {
	input {
		File sampleList
		String fastq_dir
		Boolean paired
		String bwa_index
		String chromNoScaffold
		String blacklist
		String chromosome_sizes
		String? PeakcallingControl
	}
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
		String sampleName=sample[0]
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
				BWAIndex=bwa_index,
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
		call filterTasks.remove_blacklist {
			input:
				bam=remove_scaffolds.bam_noScaffold,
				blacklist=blacklist,
				sample_name=sample_out_dir+"/"+sampleName
		}
		call filterTasks.sort_bam {
			input:
				bam=remove_blacklist.bam_noBlacklist,
				sample_name=sample_out_dir+"/"+sampleName
		}
		call filterTasks.remove_duplicates {
			input:
				bam=sort_bam.bam_sorted,
				sample_name=sample_out_dir+"/"+sampleName
		}
		call filterTasks.index_bam {
			input:
				bam=remove_duplicates.bam_noDuplicate,
				sample_name=sample_out_dir+"/"+sampleName
		}
		call pcTasks.macs2 {
			input:
				sampleName=sampleName,
				sample_out_dir=sample_out_dir,
				bam=remove_duplicates.bam_noDuplicate,
				control_bam=PeakcallingControl	
		}
		call bwTasks.read_count {
			input:
				bam=remove_duplicates.bam_noDuplicate
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
				sample_name=sample_out_dir+"/"+sampleName
		}
		call bwTasks.bedgraph_to_bigwig {
			input:
				bedgraph=bam_to_bedgraph.bedgraph,
				chromosome_sizes=chromosome_sizes,
				sample_name=sample_out_dir+"/"+sampleName
		}
	}
	output {
		Array[String] finalBam = sort_bam.bam_sorted
		Array[String] finalBamIndex = index_bam.bam_index
		Array[String] peaks = macs2.narrowPeak
		Array[String] bw = bedgraph_to_bigwig.bw
	}
}
