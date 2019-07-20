# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.9")
library("dada2")
PROJECT_DIR = "/home/rpetit/projects/supernova-16s"
fastq_path = paste0(PROJECT_DIR, "/data/fastq")
fastq_r1 <- sort(list.files(fastq_path, pattern="_R1_001.fastq", full.names = TRUE))
fastq_r2 <- sort(list.files(fastq_path, pattern="_R2_001.fastq", full.names = TRUE))
sample_names <- sapply(strsplit(basename(fastq_r1), "_"), `[`, 1)

# Quality plots
plotQualityProfile(fastq_r1[1:2])
plotQualityProfile(fastq_r2[1:2])

# Filter and Trim
filt_r1 <- file.path(PROJECT_DIR, "data", "filtered", paste0(sample_names, "_R1_filt.fastq.gz"))
filt_r2 <- file.path(PROJECT_DIR, "data", "filtered", paste0(sample_names, "_R2_filt.fastq.gz"))
names(filt_r1) <- sample_names
names(filt_r2) <- sample_names
out <- filterAndTrim(fastq_r1, filt_r1, fastq_r2, filt_r2, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learn Error Rates
error_r1 <- learnErrors(filt_r1, multithread=TRUE)
error_r2 <- learnErrors(filt_r2, multithread=TRUE)

plotErrors(error_r1, nominalQ=TRUE)

# Sample Inference
dada_r1 <- dada(filt_r1, err=error_r1, multithread=TRUE)
dada_r2 <- dada(filt_r2, err=error_r2, multithread=TRUE)

# Merge Reads
merged_reads <- mergePairs(dada_r1, filt_r1, dada_r2, filt_r2, verbose=TRUE)
head(merged_reads[[1]])

# Sequence Table
sequence_table <- makeSequenceTable(merged_reads)
dim(sequence_table)
table(nchar(getSequences(sequence_table)))

# Remove Chimeras
sequence_table_nochim <- removeBimeraDenovo(sequence_table, method="consensus", multithread=TRUE, verbose=TRUE)
dim(sequence_table_nochim)
sum(sequence_table_nochim)/sum(sequence_table)

# Track Reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_r1, getN), sapply(dada_r2, getN), sapply(merged_reads, getN), rowSums(sequence_table_nochim))
colnames(track) <- c("input", "filtered", "denoisedR1", "denoisedR2", "merged", "nonchim")
rownames(track) <- sample_names
head(track)

# Assign Taxonomy
taxa <- assignTaxonomy(sequence_table_nochim, paste0(PROJECT_DIR, "/data/reference/silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


                       
                       

