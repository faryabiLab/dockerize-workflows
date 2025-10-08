version 1.0

task macs2 {
  input {
    # --- Required ---
    String sampleName
    File bam
    File? control_bam

    # --- Runtime / container ---
    String Dockerhub_Pull = "faryabilab/macs2:0.1.0"
    Int cpu = 12
    Int mem = 25

    # --- Peak-calling parameters (all visible in inputs.json) ---
    Float? p_value
    Float q_value = 0.05

    Int? shift
    Int? extension_size
    Int? tag_size
    Int? cutoff_analysis
    Int? downsample

    Boolean call_summits = false
    Boolean fix_bimodal = false
    Boolean no_model = false
    Boolean no_lambda = false
    Boolean broad_peaks = false
    Boolean save_frag_pileup = false
  }

  command {
    set -euo pipefail

    macs2 callpeak \
      -t ${bam} \
      ~{if defined(control_bam) then "-c " + control_bam else ""} \
      -n ${sampleName} \
      ~{if defined(p_value) then "-p " + p_value else "-q " + q_value} \
      ~{if defined(shift) then "--shift " + shift else ""} \
      ~{if defined(extension_size) then "--extsize " + extension_size else ""} \
      ~{if defined(tag_size) then "--tag-size " + tag_size else ""} \
      ~{if call_summits then "--call-summits" else ""} \
      ~{if fix_bimodal then "--fix-bimodal" else ""} \
      ~{if no_model then "--nomodel" else ""} \
      ~{if save_frag_pileup then "--save-frag-pileup" else ""} \
      ~{if defined(downsample) then "--down-sample " + downsample else ""} \
      ~{if no_lambda then "--nolambda" else ""} \
      ~{if broad_peaks then "--broad" else ""} \
      ~{if defined(cutoff_analysis) then "--cutoff-analysis " + cutoff_analysis else ""} \
      --bw 300 \
      --seed 0 \
      --keep-dup 1
  }
  output {
    File narrowPeak = "${sampleName}_peaks.narrowPeak"
    File summits = "${sampleName}_summits.bed"
  }
  runtime {
    docker: Dockerhub_Pull
    cpu: cpu
    memory: mem + "G"
  }
}

task SEACR {
	input {
		String sampleName
                String? sample_out_dir
		String? Dockerhub_Pull = "faryabilab/seacr:0.1.0"

                String bedgraph
		String? control_bedgraph
		Float? top_peak_fraction

		String Normalization = "norm"
		String RunMode = "stringent"

		Int cpu = 12
		Int mem = 25

		String type

	}
	command {
		bash "/tmp/SEACR-1.3/SEACR_1.3.sh" \
		"${bedgraph}" \
		~{if defined(control_bedgraph) then control_bedgraph else top_peak_fraction} \
		"${Normalization}" \
		"${RunMode}" \
		"${sampleName}.${type}"
	}
	output {
		File seacr_out = "${sampleName}.${type}.${RunMode}.bed"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}

task bamToBedgraph {
	input {
		String? bam
		String sampleName
		String type
		String? Dockerhub_Pull = "faryabilab/bedtools:0.1.0"

		Int cpu = 12
		Int mem = 25
	}
	command {
		bedtools genomecov \
		-ibam "${bam}" \
		-bg \
		> "${sampleName}.seacr.bg.${type}"
	}
	output {
		File seacr_bg = "${sampleName}.seacr.bg.${type}"
	}
	runtime {
		docker: "${Dockerhub_Pull}"
		cpu: "${cpu}"
		mem: "${mem}"
	}
}
