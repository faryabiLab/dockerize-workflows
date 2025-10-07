setwd('/dobby/noah/dev/dockerize-workflows/benchmark/validation')

new.counts <- read.table('/dobby/noah/dev/dockerize-workflows/workflows/rna_seq/outputs/1M_SRR9336468.counts',sep="\t",header=T)
old.counts <- read.table('/dobby/noah/dev/old_workflows/rna_seq/processed/1M_SRR9336468/1M_SRR9336468.counts',sep="\t",header=T)
counts <- data.frame('old'=old.counts$X.dobby.noah.dev.old_workflows.rna_seq..processed.1M_SRR9336468.1M_SRR9336468.f.bam, 
                     'new'=new.counts$X.cromwell.executions.RNAseq.b72d5e4e.d0f6.49e8.9da9.c67abadef53e.call.count_reads_paired.shard.0.inputs.1371479256.1M_SRR9336468.sorted.bam)

M <- log2(counts$new + 1) - log2(counts$old + 1) # log fold change
A <- 0.5 * (log2(counts$new + 1) + log2(counts$old + 1)) # average intensity

pdf("rna_seq/RNAseq_OldvsNew_MAplot.pdf")
plot(A, M,
     pch = 16, cex = 0.4,
     xlab = "Average log2(count)",
     ylab = "log2(New / Old)",
     main = "MA Plot: New vs Old Counts")
abline(h = 0, col = "red")
dev.off()
