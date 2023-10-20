version 1.0

## import all relevant tasks here
import "wdl_tasks/trim_galore.wdl" as trimTasks
import "wdl_tasks/star.wdl" as starTasks
import "wdl_tasks/filter_steps.wdl" as filterTasks
import "wdl_tasks/feature_count.wdl" as quantifyTasks
import "wdl_tasks/make_bigWig.wdl" as bwTasks

workflow rnaseq {
        input {
		Boolean paired
                String fastq_dir
		String sample_out_dir
                String sampleName       
                String star_index
                String chromNoScaffold
		String? blacklist
		String GeneAnnotationFile
		String chromosome_sizes
        }

	call trimTasks.fastqc_trim {
		input:
			fastq_dir=fastq_dir,
			sample_out_dir=sample_out_dir,
			sampleName=sampleName,
			paired=paired
	}
	call starTasks.STAR {
		input:
			fastq_trimmed_single = fastqc_trim.out_fqc,
			fastq1_trimmed = fastqc_trim.out_fqc1,
			fastq2_trimmed = fastqc_trim.out_fqc2,
			star_index=star_index,
			sample_name=sampleName,
			paired=paired
	}
		
	call filterTasks.remove_scaffolds {
		input:
			bam=STAR.bam,
			chrom_no_scaff=chromNoScaffold,
			sample_name=sample_out_dir+"/"+sampleName
	}
	call filterTasks.remove_duplicates {
		input:
			bam=remove_scaffolds.bam_noScaffold,
			sample_name=sample_out_dir+"/"+sampleName,
	}
	String next = remove_duplicates.bam_noDuplicate
	if (defined(blacklist)) {
		call filterTasks.remove_blacklist {
			input:
				bam=remove_duplicates.bam_noDuplicate,
				blacklist=blacklist,
				sample_name=sample_out_dir+"/"+sampleName
		}
		String? next = remove_blacklist.bam_noBlacklist
	}
	call filterTasks.sort_bam {
		input:
			bam=next,
			sample_name=sample_out_dir+"/"+sampleName
	}
	call filterTasks.index_bam {
		input:
			bam=sort_bam.bam_sorted,
			sample_name=sample_out_dir+"/"+sampleName
	}
	if (paired) {
		call quantifyTasks.count_reads_paired {
			input:
				bam=sort_bam.bam_sorted,
				GeneAnnotationFile=GeneAnnotationFile,
				sample_name=sample_out_dir+"/"+sampleName
		}
	}
	if (!paired) {
                call quantifyTasks.count_reads_single {
                        input:
                                bam=sort_bam.bam_sorted,
                                GeneAnnotationFile=GeneAnnotationFile,
				sample_name=sample_out_dir+"/"+sampleName
                }
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
			sample_name=sample_out_dir+"/"+sampleName
	}
	call bwTasks.bedgraph_to_bigwig {
		input:
			bedgraph=bam_to_bedgraph.bedgraph,
			chromosome_sizes=chromosome_sizes,
			sample_name=sample_out_dir+"/"+sampleName
	}
	output {
		File finalBam = sort_bam.bam_sorted
		File finalBamIndex = index_bam.bam_index
		File? countsSingle = count_reads_single.gene_counts
		File? countsPaired = count_reads_paired.gene_counts
		File bw = bedgraph_to_bigwig.bw
	}
}
