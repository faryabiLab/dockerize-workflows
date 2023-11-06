version 1.0

import "wdl_tasks/trim_galore.wdl" as trimTasks
import "wdl_tasks/bwa.wdl" as bwaTasks
import "wdl_tasks/filter_steps.wdl" as filterTasks
import "wdl_tasks/feature_count.wdl" as quantTasks
import "wdl_tasks/peak_calling.wdl" as pcTasks
import "wdl_tasks/make_bigWig.wdl" as makeBWWorkflow

workflow cut_and_run {
	input {
		File sampleList
		String fastq_dir
		String project_out_dir
		Boolean paired 
		String BWAIndex
		String PeakCaller
		String? PeakCallingControl
		String ChromNoScaffold
		String Blacklist
		String ChromosomeSizes
		String? FastqSuffix
		Int? size_high
		Int? size_low
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
				fastq_suffix=FastqSuffix,
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
				BWAIndex=BWAIndex
		}
		call filterTasks.sam_to_bam {
			input:
				sam=BWA.rawSam,
				sample_name="${sample_out_dir}"+"/"+"${sampleName}"+".raw"
		}
		call filterTasks.remove_scaffolds {
			input:
				bam=sam_to_bam.bam,
				chrom_no_scaff=ChromNoScaffold,
				sample_name="${sample_out_dir}"+"/"+"${sampleName}"
		}
		call filterTasks.remove_blacklist {
			input:
				bam=remove_scaffolds.bam_noScaffold,
				blacklist=Blacklist,
				sample_name="${sample_out_dir}"+"/"+"${sampleName}"
		}
		call filterTasks.filter_discordant_pairs {
			input:
				bam=remove_blacklist.bam_noBlacklist,
				sample_name="${sample_out_dir}"+"/"+"${sampleName}"
		}
		call filterTasks.sort_bam {
			input:
				bam=filter_discordant_pairs.bam_pairedReads,
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
				threshold_hi=size_high,
				threshold_low=size_low
		}
		if ("${PeakCaller}" == "macs2") {
			call pcTasks.macs2 as allReadsPeakCalling_MACS2 {
				input:
					sampleName=sampleName,
					sample_out_dir=sample_out_dir,
					bam=remove_duplicates.bam_noDuplicate,
					control_bam=PeakCallingControl
			}
			if ("${size_high}" != "") {
				call pcTasks.macs2 as lowReadsPeakCalling_MACS2 {
					input:
						sampleName=sampleName,
						sample_out_dir=sample_out_dir,
						bam=size_filter_bam.hi,
						control_bam=PeakCallingControl
				}
			}
			if ("${size_low}" != "") {
				call pcTasks.macs2 as highReadsPeakCalling_MACS2 {
					input:
						sampleName=sampleName,
						sample_out_dir=sample_out_dir,
						bam=size_filter_bam.low,
						control_bam=PeakCallingControl
				}
			}
			
		}
		if ("${PeakCaller}" == "seacr") {
			call pcTasks.bamToBedgraph as allReads_BAMtoBG{
				input:
					bam=remove_duplicates.bam_noDuplicate,
					sampleName="${sample_out_dir}/${sampleName}",
					type="all"
			}	
			call pcTasks.SEACR as allReadsPeakCalling_SEACR{
				input:
					sampleName=sampleName,
					sample_out_dir=sample_out_dir,
					bedgraph=allReads_BAMtoBG.seacr_bg,
					control_bedgraph=PeakCallingControl,
					type="all"
			}
			if ("${size_high}" != "") {
				call pcTasks.bamToBedgraph as lowReads_BAMtoBG {
					input:
						bam=size_filter_bam.hi,
						sampleName="${sample_out_dir}/${sampleName}",
						type="low"
				}
				call pcTasks.SEACR as lowReadsPeakCalling_SEACR {
					input:
						sampleName=sampleName,
						sample_out_dir=sample_out_dir,
						bedgraph=lowReads_BAMtoBG.seacr_bg,
						control_bedgraph=PeakCallingControl,
						type="low"
				}
			}
			if ("${size_low}" != "") {
				call pcTasks.bamToBedgraph as highReads_BAMtoBG {
					input:
						bam=size_filter_bam.low,
						sampleName="${sample_out_dir}/${sampleName}",
						type="high"
				}
				call pcTasks.SEACR as highReadsPeakCalling_SEACR {
					input:
						sampleName=sampleName,
						sample_out_dir=sample_out_dir,
						bedgraph=highReads_BAMtoBG.seacr_bg,
						control_bedgraph=PeakCallingControl,
						type="high"
				}
			}
		}
		call makeBWWorkflow.makeBigWig as BW_allReads {
			input:
				bam=remove_duplicates.bam_noDuplicate,
				chromosome_sizes=ChromosomeSizes,
				sampleName="${sample_out_dir}"+"/"+"${sampleName}.AllReads"
		}
		if ("${size_high}" != "") {
			call makeBWWorkflow.makeBigWig as BW_hiThresh {
				input:
					bam=size_filter_bam.hi,
					chromosome_sizes=ChromosomeSizes,
					sampleName="${sample_out_dir}"+"/"+"${sampleName}.HighThreshold"
			}
		}
	}
	output {
		Array[String] bam = remove_duplicates.bam_noDuplicate
		Array[String] bw = BW_allReads.bw
		Array[String?] SeacrAllReadsPeaks = allReadsPeakCalling_SEACR.seacr_out
		Array[String?] SeacrLowReads = lowReadsPeakCalling_SEACR.seacr_out
		Array[String?] SeacrHighReads = highReadsPeakCalling_SEACR.seacr_out
		Array[String?] AllReadsPeaks = allReadsPeakCalling_MACS2.narrowPeak
		Array[String?] LowReadsPeaks = lowReadsPeakCalling_MACS2.narrowPeak
		Array[String?] HighReadsPeaks = highReadsPeakCalling_MACS2.narrowPeak
		Array[String?] bam_low = size_filter_bam.low
		Array[String?] bam_hi = size_filter_bam.hi
		Array[String?] bw_hiThresh = BW_hiThresh.bw
		#File? bw_lowThresh = BW_lowThresh.bw
	}
}
