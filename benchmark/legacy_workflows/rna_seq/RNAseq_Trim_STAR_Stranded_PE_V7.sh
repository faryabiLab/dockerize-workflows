
#### for paired-end RNAseq alignment with STAR on Plutus
#### using our old parameters
#### choice of one/two/no trimming
#### featureCount on no rRNA, no chrM gtf.

#### input ----------------------------------------------------------------------------------//
#### two adaptors
if [ $# -eq 7 ]; then
        SPECIES=${1}
        STARIND=${2}
        SEQDIR=${3}
        OUTDIR=${4}
        EXTENSION=${5}
        ADAPT_01=${6}
        ADAPT_02=${7}
#### one adaptor
elif [ $# -eq 6 ]; then
        SPECIES=${1} 
        STARIND=${2}
        SEQDIR=${3}
        OUTDIR=${4}
        EXTENSION=${5}
        ADAPT_01=${6}
#### no adaptor
elif [ $# -eq 5 ]; then
        SPECIES=${1}  
        STARIND=${2}
        SEQDIR=${3}
        OUTDIR=${4}
        EXTENSION=${5}
else     
        echo -e "This is the script for STAR RNA-seq alignment. \n 
        Usage: 
                ./RNAseq_Trim_STAR_Stranded_PE_V7.sh <SPECIES=hg19/mm10> <STARIndexDIR> <SequenceDIR> <AlignDIR> <SequenceExtension> <OPT1: adaptor1> <OPT2: adaptor2> \n
        "     
        exit 1
fi
  
#### choose reference genome and annotation -------------------------------------------------//
if [ ${SPECIES} == "hg19" ]; then
  
        CHROM="/mnt/data0/fdb/igenome/Homo_sapiens/Ensembl/GRCh37.75/Annotation/Genes"
        CHROMGTF="${CHROM}/Homo_sapiens_GRCh37_75_noRRNA_noM.gtf"
        CHROMINFO="${CHROM}/ChromInfo.txt"
        CHROMCHR="${CHROM}/ChromNoPatch.bed"
  
elif [ ${SPECIES} == "mm10" ]; then
  
        CHROM="/mnt/data0/fdb/igenome/Mus_musculus/Ensembl/GRCm38.87/Annotation/Genes"
        CHROMGTF="${CHROM}/Mus_musculus.GRCm38.87.noRRNA.noM.gtf"
        CHROMINFO="${CHROM}/ChromInfo.txt"
        CHROMCHR="${CHROM}/ChromNoPatch.bed"
elif [ ${SPECIES} == "yeast" ]; then
        CHROM="/dobby/noah/dev/dockerize-workflows/test_data"
        CHROMGTF="${CHROM}/Saccharomyces_cerevisiae.R64-1-1.gtf"
        CHROMINFO="${CHROM}/ChromSizes.txt"
        CHROMCHR="${CHROM}/ChromNoPatch.txt"
else     
        echo "unrecognized species"
        exit 1
fi
                       
#### filter:  read unmapped; read fails platform/vendor quality checks
CPU=6
FILTER=516
MEM=15000000000

trim_glore_ext1="val_1.fq.gz"
trim_glore_ext2="val_2.fq.gz"

########################################################################
for FILE in ${SEQDIR}/*${EXTENSION}
############### loop begins #######
do               
                 
#### file name without path
FORWARD=${FILE##*/}
#### find the second end of the pair
EXTENSION_02=$(echo ${EXTENSION} | sed -e "s/1/2/")
BACKWARD=`echo ${FORWARD} | sed -e "s/${EXTENSION}/${EXTENSION_02}/g"`
                 
#p1=`echo $(basename "${EXTENSION}" ".fastq.gz")`
#p2=`echo $(basename "${EXTENSION_02}" ".fastq.gz")`
                 
#### prefix shared by two ends
#kb=$(basename "$FORWARD" _1.fastq.gz)
kb=`echo $(basename "$FORWARD" ${EXTENSION})`
                 
#### trim galore output 
                 
#### shell script name
INPUTFQ=${kb}    
                 
#### trim galore output name
TRIM_01=$(echo ${FORWARD} | cut -f 1 -d ".")_${trim_glore_ext1}
TRIM_02=$(echo ${BACKWARD} | cut -f 1 -d ".")_${trim_glore_ext2}
                 
#### full working path
PREFIX=${OUTDIR}/processed/${kb}
mkdir -p ${PREFIX}

#### output extensions 
OUTBAMRAW=${kb}.raw.bam
OUTBAM=${kb}.bam
OUTBAMDM_NS=${kb}_unsort.dm.bam
OUTBAMDM=${kb}.dm.bam
OUTBAMDM_METRIC=${kb}.dupmetric.bam
COUNT=${kb}.counts

#### final output
OUTBAM_FILTER=${kb}.f.bam
BGNNJ=${kb}_NJ_norm_unsort.bg
BWNNJ=${kb}_NJ_norm.bw

echo "
#!/bin/bash" >> ${INPUTFQ}.sh

echo " cd ${PREFIX} " >> ${INPUTFQ}.sh

#### trim adapter and define the input to STAR ---------------------------------------------------------------------// 
if [ ! -z ${ADAPT_02} ] && [ ! -z ${ADAPT_01} ];then  
        #### trim
        echo "two adaptors"  >> ${INPUTFQ}.sh
        echo "trim_galore --fastqc -q 15 --phred33 --gzip --stringency 5 -e 0.1 --length 20 -a ${ADAPT_01} -a2 ${ADAPT_02} --paired ${SEQDIR}/${FORWARD} ${SEQDIR}/${BACKWARD} -o ${PREFIX}" >> ${INPUTFQ}.sh
elif [ ! -z ${ADAPT_01} ];then
        #### trim
        echo "one adaptor"
        echo "trim_galore --fastqc -q 15 --phred33 --gzip --stringency 5 -e 0.1 --length 20 -a ${ADAPT_01} --paired ${SEQDIR}/${FORWARD} ${SEQDIR}/${BACKWARD} -o ${PREFIX}" >> ${INPUTFQ}.sh
elif [ -d "${SEQDIR}" ] && [ -d "${OUTDIR}" ] && [ -n ${EXTENSION} ]; then
        printf "${FORWARD} and ${BACKWARD} \n" 
        echo "default trimming"
        echo "trim_galore --fastqc -q 15 --phred33 --gzip --stringency 5 -e 0.1 --length 20 --illumina --paired ${SEQDIR}/${FORWARD} ${SEQDIR}/${BACKWARD} -o ${PREFIX}" >> ${INPUTFQ}.sh
fi

#### STAR -----------------------------------------------------------------------------------------------------------------------------------------------//
echo " ## run the star and make the aligned bam ## " >> ${INPUTFQ}.sh

echo "echo \"######## beginning alignment ########\"" >> ${INPUTFQ}.sh

echo " STAR \\
--genomeDir ${STARIND} \\
--outFilterType BySJout \\
--readFilesIn ${PREFIX}/${TRIM_01} ${PREFIX}/${TRIM_02} \\
--readFilesCommand zcat \\
--runThreadN ${CPU} \\
--outSAMattributes Standard \\
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \\
--alignIntronMax 100000 \\
--outSAMstrandField intronMotif \\
--outSAMattrRGline ID:1 SM:${kb} LB:lib1 PL:ILLUMINA PU:unit1 \\
--outFileNamePrefix ${PREFIX}/ \\
--outSAMunmapped Within \\
--chimSegmentMin 25 \\
--chimJunctionOverhangMin 25 \\
--outStd SAM | samtools view -bS - | samtools sort -m ${MEM} - -o ${PREFIX}/${OUTBAMRAW}" >> ${INPUTFQ}.sh

echo " " >> ${INPUTFQ}.sh
echo " ## index the aligned bam" >> ${INPUTFQ}.sh
echo "bamtools index -in ${PREFIX}/${OUTBAMRAW}" >> ${INPUTFQ}.sh

#### sort raw alignment --> remove scaffolds --> index -------------------------------------------------------------------//
echo " " >> ${INPUTFQ}.sh
echo "## filter out scaffold and index the final bam##" >> ${INPUTFQ}.sh
echo "samtools view -h -L ${CHROMCHR} ${PREFIX}/${OUTBAMRAW} | samtools sort -m ${MEM} - -o ${PREFIX}/${OUTBAM}" >> ${INPUTFQ}.sh
echo "bamtools index -in ${PREFIX}/${OUTBAM}" >> ${INPUTFQ}.sh

#### mark duplicates ----------------------------------------------------------------------------------------------------//
echo " " >> ${INPUTFQ}.sh
echo "## mark the duplicates ##" >> ${INPUTFQ}.sh
echo "picard MarkDuplicates \\
AS=true \\
M=${PREFIX}/${OUTBAMDM_METRIC} O=${PREFIX}/${OUTBAMDM_NS} I=${PREFIX}/${OUTBAM} \\
REMOVE_DUPLICATES=false \\
VALIDATION_STRINGENCY=SILENT" >> ${INPUTFQ}.sh

echo " " >> ${INPUTFQ}.sh
echo "## index the duplication marked bam file" >> ${INPUTFQ}.sh
echo "samtools sort -m ${MEM} ${PREFIX}/${OUTBAMDM_NS} -o ${PREFIX}/${OUTBAMDM}" >> ${INPUTFQ}.sh
echo "bamtools index -in ${PREFIX}/${OUTBAMDM}" >> ${INPUTFQ}.sh

#### filter -------------------------------------------------------------------------------------------------------------// 
echo "# filter by read unmapped, fails platform quality, PCR duplicates" >> ${INPUTFQ}.sh
echo "samtools view -b -F ${FILTER} ${PREFIX}/${OUTBAMDM} \\
| samtools sort -m ${MEM} - -o ${PREFIX}/${OUTBAM_FILTER}" >> ${INPUTFQ}.sh

echo " " >> ${INPUTFQ}.sh
echo "bamtools index -in ${PREFIX}/${OUTBAM_FILTER}" >> ${INPUTFQ}.sh

#### get the transcripts counts -----------------------------------------------------------------------------------------//
echo " " >> ${INPUTFQ}.sh
echo "## rna count ##" >> ${INPUTFQ}.sh
echo "echo \"######## read count ########\"" >> ${INPUTFQ}.sh

echo "featureCounts -p -t exon -g gene_id -s 1 -O -T ${CPU} -a ${CHROMGTF} -o ${PREFIX}/${COUNT} ${PREFIX}/${OUTBAM_FILTER}" >> ${INPUTFQ}.sh
echo " " >> ${INPUTFQ}.sh

#### get normalized BigWig ----------------------------------------------------------------------------------------------//
echo " " >> ${INPUTFQ}.sh   
echo "echo \"######## beginning normalized bigWigg files ########\"" >> ${INPUTFQ}.sh
echo "#### normalizing ##" >> ${INPUTFQ}.sh
echo "T=\$(samtools view -c ${PREFIX}/${OUTBAM_FILTER})" >> ${INPUTFQ}.sh
echo "FACTOR=\`echo \"scale=10; 1000000 / \${T}\" | bc -l\`" >> ${INPUTFQ}.sh
 
#### Bed12 is necessary for making correct junctions
#### Unfortunately we can't pipe genomecoverage directly into bedGraphToBigWig because bedGraphToBigWig does a few passes over the file
echo "bamToBed -i ${PREFIX}/${OUTBAM_FILTER} -bed12 | bed12ToBed6 -i stdin | genomeCoverageBed -bg -i stdin -g ${CHROMINFO} -scale \${FACTOR} | sort -k1,1 -k2,2n > ${PREFIX}/${BGNNJ}" >> ${INPUTFQ}.sh  
echo "bedGraphToBigWig ${PREFIX}/${BGNNJ} ${CHROMINFO} ${PREFIX}/${BWNNJ}" >> ${INPUTFQ}.sh  
 
chmod 755 ${INPUTFQ}.sh

echo "rm *bg" >> ${INPUTFQ}.sh
echo "rm tmp.bam" >> ${INPUTFQ}.sh
 
done
