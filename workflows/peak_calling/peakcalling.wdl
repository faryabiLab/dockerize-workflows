version 1.0

import "wdl_tasks/peak_calling.wdl" as peakcall

workflow MACS2_PeakCalling {
    meta {
        description: "Boilerplate workflow for MACS2 peak calling"
    }

    input {
        # Required inputs
        String sampleList
        Boolean paired_end = false

        # Optional control BAM
        File? control_bam

        # Parameters
        String genome_size = "hs"  # hs=human, mm=mouse, ce=worm, dm=fly
        Float q_value = 0.05
        Boolean call_summits = true
        String peak_type = "narrow"  # or "broad"

        # Runtime / Docker
        Int cpu = 4
        Int mem = 8
    }
    
    Array[Array[String]] samples = read_tsv(sampleList)
    scatter (sample in samples){

        String sample_id = sample[0]	
        String bam = sample[1]     

        call peakcall.MACS2_CallPeaks {
            input:
		treatment_bam = bam,
		sample_name = sample_id,
                control_bam = control_bam,
                genome_size = genome_size,
                q_value = q_value,
                call_summits = call_summits,
                paired_end = paired_end,
                peak_type = peak_type,
                cpu = cpu,
                mem = mem
        }
    }
    output {
        Array[File] narrowPeak = MACS2_CallPeaks.narrowPeak
        Array[File?] summits = MACS2_CallPeaks.summits
        Array[File] xls = MACS2_CallPeaks.xls
    }
}
