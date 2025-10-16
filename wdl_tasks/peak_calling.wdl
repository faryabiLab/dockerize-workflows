version 1.0

task MACS2_CallPeaks 
{
    input 
    {
        String sample_name
        File treatment_bam
        File? control_bam
        String genome_size
        Float q_value
        Boolean call_summits
        Boolean paired_end
        String peak_type
        String DockerhubPull = "faryabilab/macs2:0.1.0"
        Int cpu = 4
        Int mem = 8
     }

    command <<<
        set -euo pipefail
        mkdir -p macs2_out

        # Build optional flags safely
        FLAGS=""
        if [[ "~{paired_end}" == "true" ]]; then FLAGS="$FLAGS --format BAMPE"; else FLAGS="$FLAGS --format BAM"; fi
        if [[ "~{call_summits}" == "true" ]]; then FLAGS="$FLAGS --call-summits"; fi
        if [[ "~{peak_type}" == "broad" ]]; then FLAGS="$FLAGS --broad"; fi
        CONTROL_FLAG=""
        if [[ -n "~{control_bam}" && "~{control_bam}" != "null" ]]; then CONTROL_FLAG="-c ~{control_bam}"; fi

        macs2 callpeak \
            -t ~{treatment_bam} \
            $CONTROL_FLAG \
            -n ~{sample_name} \
            -g ~{genome_size} \
            --qvalue ~{q_value} \
            --outdir macs2_out \
            $FLAGS

        echo "MACS2 version: $(macs2 --version)" > macs2_out/run.log
        echo "Genome size: ~{genome_size}" >> macs2_out/run.log
        echo "Q-value cutoff: ~{q_value}" >> macs2_out/run.log
    >>>

    output {
        File narrowPeak = "${sample_name}_peaks.narrowPeak"
        File? summits = "${sample_name}_summits.bed"
        File xls = "${sample_name}_peaks.xls"
    }

    runtime {
        docker: "${DockerhubPull}"
        cpu: "${cpu}"
        mem: "${mem}"
    }
}

