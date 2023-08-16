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
                File fastq1
                File? fastq2

                String sampleName
                # trim_galore args
                String trim_out_dir
                Int qualityCutoff
                Int stringencyCutoff
                Float errorRate
                Int lengthCutoff
                # STAR args             
                String star_index
                String star_out_dir
                String outFilterType
                String readFilesCommand 
                String outSamAttributes 
                String outFilterIntronMotifs 
                Int alignIntronMax
                String outSamstrandField 
                String outSAMunmapped
                Int chimSegmentMin 
                Int chimJunctionOverhangMin 
                String outSAMtype 
                # Remove scaffold args
                File chromNoScaffold
		String removeScaffold_out_dir
		# Remove duplicates args
		String removeDuplicate_out_dir
		String PicardRemoveDuplicates
		String PicardValidationStringency
		String PicardMetricsFile
		# Remove blacklist args
		String removeBlacklist_out_dir
		File blacklist
		# Sort bam args
		String sortBam_out_dir
		# Index bam args
		String indexBam_out_dir
		# Quantify args
		String counts_out_dir
		File GeneAnnotationFile
		String AttributeType
		String GTFAttributeType
		String Stranded
		# Make BigWig args
		String bw_out_dir
		File chromosome_sizes
        }
        call trimTasks.fastqc_trim {
                input:
                        fastq1=fastq1,
                        fastq2=fastq2,
                        sampleName=sampleName,
                        quality=qualityCutoff,
                        stringency=stringencyCutoff,
                        e=errorRate,
                        length=lengthCutoff
        }
        call starTasks.STAR {
                input:
                        fastq1_trimmed=fastqc_trim.out_fqc1,
                        fastq2_trimmed=fastqc_trim.out_fqc2,
                        star_index=star_index,
                        out_dir=out_dir,
                        sample_name=sampleName,
                        readFilesCommand=readFilesCommand,
                        outSamAttrivutes=outSamAttributes,
                        outFilterIntronMotifs=outFilterIntronMotifs,
			alignIntronMax=alignIntronMax,
			outSamstrandField=outSamstrandField,
			outSAMunmapped=outSAMunmapped,
			chimSegmentMin=chimSegmentMin,
			ChimJunctionOverhangMin=chimJunctionOverhangMin,
			outSAMtype=outSAMtype
	}
	call filterTasks.remove_scaffolds {
		input:
			bam=STAR.bam,
			chrom_no_scaff=chromNoScaffold,
			out_dir=removeScaffold_out_dir,
			sample_name=sampleName
	}
	call filterTasks.remove_duplicates {
		input:
			bam=remove_scaffolds.bam_sorted,
			out_dir=removeDuplicate_out_dir,
			sample_name=sampleName,
			PicardRemoveDuplicates=PicardRemoveDuplicates,
			PicardValidationStringency=PicardValidationStringency,
			PicardMetricsFile=PicardMetricsFile
	}
	call filterTasks.remove_blacklist {
		input:
			bam=remove_duplicates.bam_noscaff,
			blacklist=blacklist,
			out_dir=removeBlacklist_out_dir,
			sample_name=sampleName
	}
	call filterTasks.sort_bam {
		input:
			bam=remove_blacklist.bam_noBlacklist,
			out_dir=sortBam_out_dir,
			sample_name=sampleName
	}
	call filterTasks.index_bam {
		input:
			bam=sort_bam.bam_sorted,
			out_dir=indexBam_out_dir,
			sample_name=sampleName
	}
	if (paired) {
		call quantifyTasks.count_reads_paired {
			input:
				bam=sort_bam.bam_sorted,
				GeneAnnotationFile=GeneAnnotationFile,
				GTFAttributeType=GTFAttributeType,
				Stranded=Stranded,
				out_dir=counts_out_dir
		}
	}
	if (!paired) {
                call quantifyTasks.count_reads_single {
                        input:
                                bam=sort_bam.bam_sorted,
                                GeneAnnotationFile=GeneAnnotationFile,
                                GTFAttributeType=GTFAttributeType,
                                Stranded=Stranded,
                                out_dir=counts_out_dir
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
			out_dir=bw_out_dir,
			sample_name=sampleName
	}
	call bwTasks.bedgraph_to_bigwig {
		input:
			bedgraph=bam_to_bedgraph.bedgraph,
			chromosome_sizes=chromosome_sizes,
			out_dir=bw_out_dir,
			sample_name=sampleName
	}
	output {
		File finalBam = sort_bam.bam_sorted
		File finalBamIndex = index_bam.bam_index
		File? countsSingle = count_reads_single.gene_counts
		File? countsPaired = count_reads_paired.gene_counts
		File bw = bedgraph_to_bigwig.bw
	}
}



