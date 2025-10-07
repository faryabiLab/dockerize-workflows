#!/usr/bin/env bash

#### 17/01/10
#### Yeqiao

#### script for creating STAR alignment file for each stranded paired-end sample
#### WARNING: please name all paired samples as *_R1.fastq.gz, #_R2.fastq.gz

#### FIXED FOR PLUTUS ------------------------------------------------------------
#### script for alignment 
#STARSCRIPT=/mnt/data0/yeqiao/analysis/pipeline/RNA_align_STAR/170111_STAR_PE/RNAseq_Trim_STAR_PE_v6.sh
STARSCRIPT=/dobby/noah/dev/old_workflows/rna_seq/RNAseq_Trim_STAR_Stranded_PE_V7.sh

#### CHANGE FOR EACH EXPERIMENT ---------------------------------------------------
#### species: hg19 or mm10
SPECIES=yeast
 
#### STAR index 
STARIND="/dobby/noah/dev/dockerize-workflows/test_data/STAR_idx"

#### directory to fastq files
SEQDIR="/dobby/noah/dev/dockerize-workflows/test_data/RNA-seq/RNA-Seq_Sample_Files" 

#### directory to output alignment
OUTDIR="/dobby/noah/dev/old_workflows/rna_seq/"

#### zipped fastq file extension
EXTENSION="_1.fastq.gz"

#### one adaptor; comment out if none.
#ADAPT_01="AGATCGGAAGAGC"

#### two adaptors; comment out if none.
#ADAPT_02="GCTCTTCCGATCT"

#### RUN SCRIPT -------------------------------------------------------------------
bash ${STARSCRIPT} ${SPECIES} ${STARIND} ${SEQDIR} ${OUTDIR} ${EXTENSION}
#bash ${STARSCRIPT} ${SEQDIR} ${OUTDIR} ${EXTENSION} ${ADAPT_01}

