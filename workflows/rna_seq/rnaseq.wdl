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

workflow RNAseq {
	input {
		File sampleList
		String fastq_dir
		String project_out_dir
		Boolean paired       
                String star_index
                String chromNoScaffold
		String blacklist
		String GeneAnnotationFile
		String chromosome_sizes
		# trim_galore
		Int quality = 15
		Int stringency = 5
		Float e = 0.1
		Int length = 20
		# STAR
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
		# remove duplicates
		String PicardRemoveDuplicates = "true"
		String PicardValidationStringency = "SILENT"
		String PicardMetricsFile = "removeDuplicate_metrics.txt"
	}	
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
		String sampleName=sample[0]
		String sample_out_dir=project_out_dir+"/"+sampleName
		call trimTasks.fastqc_trim {
			input:
				fastq_dir=fastq_dir,
				sample_out_dir=sample_out_dir,
				sampleName=sampleName,
				paired=paired,
				quality=quality,
				stringency=stringency,
				e=e,
				length=length
		}
		call starTasks.STAR {
			input:
				fastq_trimmed_single = fastqc_trim.out_fqc,
				fastq1_trimmed = fastqc_trim.out_fqc1,
				fastq2_trimmed = fastqc_trim.out_fqc2,
				star_index=star_index,
				sample_name=sample_out_dir + "/" + sampleName,
				paired=paired,
				outFilterType=outFilterType,
				readFilesCommand=readFilesCommand,
				outSamAttributes=outSamAttributes,
				outFilterIntronMotifs=outFilterIntronMotifs,
				alignIntronMax=alignIntronMax,
				outSamstrandField=outSamstrandField,
				outSAMunmapped=outSAMunmapped,
				chimSegmentMin=chimSegmentMin,
				chimJunctionOverhangMin=chimJunctionOverhangMin,
				outSAMtype=outSAMtype,
				STAR_cpu=STAR_cpu,
				STAR_mem=STAR_mem
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
				PicardRemoveDuplicates=PicardRemoveDuplicates,
				PicardValidationStringency=PicardValidationStringency,
				PicardMetricsFile=PicardMetricsFile
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
        }	
	output {
		Array[String] finalBam = sort_bam.bam_sorted
		Array[String] finalBamIndex = index_bam.bam_index
		Array[String?] countsSingle = count_reads_single.gene_counts
		Array[String?] countsPaired = count_reads_paired.gene_counts
		Array[String] bw = bedgraph_to_bigwig.bw
	}
}
