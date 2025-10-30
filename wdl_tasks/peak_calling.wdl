version 1.0

###############################################
# 1. Fragment size estimation task (optional)
###############################################
task EstimateFragSize {
  input {
    File bam
    String sample_name

    String Dockerhub_Pull = "faryabilab/homer:v1"
    Int cpu = 12
    Int mem = 50
  }

  command <<< 
    bash <<'EOF'
    FRAGSIZE=$(homer findPeaks -style factor -size given -input ~{bam} 2>/dev/null \
        | grep 'fragment length' \
        | awk '{print $NF}')
    echo "Fragment size: $FRAGSIZE"
    EOF
  >>>

  output {
    Int estimated_fragment_size = read_int("fragsize.txt")
  }

  runtime {
    docker: "${Dockerhub_Pull}"
    cpu: "${cpu}"
    mem: "${mem}"
  }
}

###############################################
# 2. MACS2 peak calling task
###############################################
task MACS2_CallPeaks {
  input {
    String Dockerhub_Pull = "faryabilab/macs2:0.1.0"

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

  String broad_flag = if (peak_type == "broad") then "--broad" else ""

  command {
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
      ~{broad_flag} \
      ~{if paired_end then "--format BAMPE" else "--format BAM"} \
      --shift ~{shift} \
      --extsize $FRAGSIZE
  }

  output {
    File narrowPeak = "${sample_name}_peaks.narrowPeak"
    File? summits = "${sample_name}_summits.bed"
    File xls = "${sample_name}_peaks.xls"
  }

  runtime {
    docker: "${Dockerhub_Pull}"
    cpu: "${cpu}"
    mem: "${mem}"
  }
}

