version 1.0

workflow MACS2_PeakCalling {
    meta {
        description: "Boilerplate workflow for MACS2 peak calling"
    }

    input {
        # Required inputs
        String sample_name
        File treatment_bam
        File treatment_bam_index
        Boolean paired_end = false

        # Optional control BAM
        File? control_bam
        File? control_bam_index

        # Parameters
        String genome_size = "hs"  # hs=human, mm=mouse, ce=worm, dm=fly
        Float q_value = 0.05
        Boolean call_summits = true
        String peak_type = "narrow"  # or "broad"

        # Runtime / Docker
        String docker_image = "faryabilab/macs2:2.2.9.1"
        Int cpu = 4
        Int memory_gb = 8
    }

    call MACS2_CallPeaks {
        input:
            sample_name = sample_name,
            treatment_bam = treatment_bam,
            control_bam = control_bam,
            genome_size = genome_size,
            q_value = q_value,
            call_summits = call_summits,
            paired_end = paired_end,
            peak_type = peak_type,
            docker_image = docker_image,
            cpu = cpu,
            memory_gb = memory_gb
    }

    output {
        File narrowPeak = MACS2_CallPeaks.narrowPeak
        File summits = MACS2_CallPeaks.summits
        File xls = MACS2_CallPeaks.xls
        File log = MACS2_CallPeaks.log
    }
}
