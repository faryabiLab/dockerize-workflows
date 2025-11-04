## Workflows
### Quality control, alignment, filtering, and track creation
Each assay-named workflow does the following general steps: \
1. Adapter trimming and Quality Control checks with `FastQC`.
2. Alignment to a reference genome (`STAR` for RNA, `BWA mem` for ChIP/ATAC)
3. Multiple rounds of filtering after alignment:
    * Remove scaffolds (keep list of cannonical chromosomes in `chromNoScaffold` input parameter
    * Mark and remove duplicate alignments
    * Remove unmapped alignments
    * Conditionally remove blacklisted regions, if such a file is supplied to the `blackList` input parameter
    * Sort and index the final bam.
4. Produce normalized (TPM) track (bigWig) files

`README`s in each of the subdirectories outline parameters, steps, and outputs of each workflow.

