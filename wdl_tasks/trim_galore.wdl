version 1.0

task fastqc_trim {
        input {
		#### REQUIRED
                String fastq_dir
		String sample_out_dir
		String sampleName
		####

		Int quality = 15
		Int stringency = 5
		Float e = 0.1
		Int length = 20
        }
	String fastq1 = fastq_dir+"/"+sampleName+"_R1.fastq.gz"
	String fastq2 = fastq_dir+"/"+sampleName+"_R2.fastq.gz"
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
                "${fastq1}" "${fastq2}"
        }
        output {
                File out_fqc1 = "${sample_out_dir}/${sampleName}"+"_val_1.fq.gz"
                File? out_fqc2 = "${sample_out_dir}/${sampleName}"+"_val_2.fq.gz"
        }
        runtime {
                docker: "faryabilab/trim_galore:0.10"
        }

}
