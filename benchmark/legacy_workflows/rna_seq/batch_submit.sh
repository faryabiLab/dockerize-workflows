#!/bin/bash

#### path to scripts
OUTDIR=/dobby/noah/dev/old_workflows/rna_seq

#### batch nohup
for FILE in ${OUTDIR}/*SRR*.sh
do
    echo "processing $FILE"
    chmod 755 ${FILE}
    #nohup sh ${FILE} > ${FILE}.out &

    # run the alignment_scripts.txt with parallel as
    #cat alignment_scripts.txt | nohup parallel -j 5 &
    echo "bash ${FILE} 2>&1 | tee ${FILE}.out" >> alignment_scripts.txt
    cat alignment_scripts.txt | nohup parallel -j 5 &
done
