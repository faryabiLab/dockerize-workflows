version 1.0

task fastqc_trim {
        input {
                File fastq1
                File? fastq2
		String sampleName
		String out_dir
		Int quality
		Int stringency
		Float e
		Int length
        }
        command {
                trim_galore -q ${quality} \
                --paired \
                --fastqc \
                --phred33 \
                --gzip \
                --basename ${sampleName} \
                --stringency ${stringency} \
                -e ${e} \
                --length ${length} \
                -o ${out_dir} \
                "${fastq1}" \
                "${fastq2}"
        }
        output {
                File out_fqc1 = "${out_dir}"+"/"+"${sampleName}"+"_val_1.fq.gz"
                File? out_fqc2 = "${out_dir}"+"/"+"${sampleName}"+"_val_2.fq.gz"
        }
        runtime {
                cpu: 10
                docker: "faryabilab/trim_galore:0.10"
        }

}
