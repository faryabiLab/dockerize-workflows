version 1.0

task fastqc {
        input {
                File fastq1
                File fastq2
                String fq_suffix
		Int quality
		Int stringency
		Float e
		Int length
        }
        String fqc = basename(fastq1, "${fq_suffix}")
        String fqc_full = basename(fastq1)
        command {
                trim_galore -q ${quality} \
                --paired \
                --fastqc \
                --phred33 \
                --gzip \
                --basename ${fqc} \
                --stringency ${stringency} \
                -e ${e} \
                --length ${length} \
                -o 01.fastqc \
                "${fastq1}" \
                "${fastq2}"
        }
        output {
                File out_fqc1 = "01.fastqc/"+"${fqc}"+"_val_1.fq.gz"
                File out_fqc2 = "01.fastqc/"+"${fqc}"+"_val_2.fq.gz"
        }
        runtime {
                cpu: 10
                docker: "faryabilab/trim_galore:0.10"
        }

}

