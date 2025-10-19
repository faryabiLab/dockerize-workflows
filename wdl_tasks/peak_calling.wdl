version 1.0

###############################################
# 1. Fragment size estimation task (optional)
###############################################
task EstimateFragSize {
  input {
    File bam
    String sample_name
  }

  command <<<
    set -euo pipefail
    echo "Estimating fragment size for ~{sample_name}..."

    # Try to estimate fragment size using HOMER or any other tool
    FRAGSIZE=$(homer findPeaks -style factor -size given -input ~{bam} 2>/dev/null | grep 'fragment length' | awk '{print $NF}')

    if [ -z "$FRAGSIZE" ]; then
      echo "WARNING: Could not estimate fragment size, defaulting to 200 bp."
      FRAGSIZE=200
    fi

    echo $FRAGSIZE > fragsize.txt
  >>>

  output {
    Int estimated_fragment_size = read_int("fragsize.txt")
  }

  runtime {
    docker: "faryabilab/homer:0.1.0"
    cpu: 2
    memory: "4G"
  }
}

###############################################
# 2. MACS2 peak calling task
###############################################
task MACS2_CallPeaks {
  input {
    # Required
    File treatment_bam
    String sample_name

    # Optional control
    File? control_bam

    # Significance thresholds (mutually exclusive)
    Float q_value = 0.05
    Float? p_value

    # Genome & peak calling parameters
    String genome_size = "hs"  # hs=human, mm=mouse, ce=worm, dm=fly
    Boolean call_summits = true
    Boolean paired_end = false
    String peak_type = "narrow"  # or "broad"

    # Fragment size handling
    Int? estimated_fragment_size
    Int default_extsize = 200
    Int shift = 0

    # Runtime
    Int cpu = 4
    Int mem = 8
  }

  command <<<
    set -euo pipefail

    # Choose fragment size (estimated if available, else default)
    FRAGSIZE=~{select_first([estimated_fragment_size, default_extsize])}
    echo "Using fragment size: $FRAGSIZE bp"

    echo "Running MACS2 peak calling..."
    macs2 callpeak \
      -t ~{treatment_bam} \
      ~{if defined(control_bam) then ("-c " + control_bam) else ""} \
      -n ~{sample_name} \
      -g ~{genome_size} \
      ~{if defined(p_value) then ("--pvalue " + p_value) else ("--qvalue " + q_value)} \
      ~{if call_summits then "--call-summits" else ""} \
      ~{if peak_type == "broad" then "--broad" else ""} \
      ~{if paired_end then "--format BAMPE" else "--format BAM"} \
      --shift ~{shift} \
      --extsize $FRAGSIZE
  >>>

  output {
    File narrowPeak = "${sample_name}_peaks.narrowPeak"
    File? summits = "${sample_name}_summits.bed"
    File xls = "${sample_name}_peaks.xls"
  }

  runtime {
    docker: "faryabilab/macs2:0.1.0"
    cpu: cpu
    memory: "~{mem}G"
  }
}

