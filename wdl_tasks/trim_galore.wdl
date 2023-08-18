version 1.0

task fastqc_trim {
        input {
                String fastq1
                String? fastq2
		String sample_out_dir
		String sampleName
		Int quality
		Int stringency
		Float e
		Int length
        }
        command { 
                trim_galore \
                -q ${quality} \
                --paired \
                --fastqc \
                --phred33 \
                --gzip \
                --stringency ${stringency} \
                -e ${e} \
                --basename ${sampleName} \
                --length ${length} \
                -o ${sample_out_dir} \
                "${fastq1}" \
                "${fastq2}"
        }
        output {
                File out_fqc1 = "${sample_out_dir}/${sampleName}"+"_val_1.fq.gz"
                File? out_fqc2 = "${sample_out_dir}/${sampleName}"+"_val_2.fq.gz"
        }
        runtime {
                docker: "faryabilab/trim_galore:0.10"
        }

}
