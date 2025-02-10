# This pipeline is for when we have unrealistically large datasets that:
#        -A) Running anything in pooled mode is unfeasible even on our HPC (and more importantly our research question does not REQUIRE this for a fine-scale ability to assess uncertainty in richness estimates)
#        -B) Processes must be run parralel

####### Load required packages #######

suppressMessages(library(dada2)); print(paste0("loading dada 2 version ", packageVersion("dada2")))
suppressMessages(library(ShortRead)); print(paste0("loading ShortRead version ", packageVersion("ShortRead")))
suppressMessages(library(Biostrings)); print(paste0("loading Biostrings 2 version ", packageVersion("Biostrings")))
suppressMessages(library(stringr)); print(paste0("loading stringr version ", packageVersion("stringr")))

####### File parsing #######

# set working directory:
wdPath <- getwd()
list.files(wdPath)

print("Defining dummy variable inputs...")

# Define the dataset dummy variables
dataset.value <- "dataset.value.input"
print(paste0("dataset.value = ", dataset.value))
run.value <- "run.value.input"
print(paste0("run.value = ", run.value))
rc.value <- "rc.value.input"
print(paste0("rc.value = ", rc.value))
concat.value <- "concat.value.input"
print(paste0("concat.value = ", concat.value))
pool.value <- "pool.value.input"
print(paste0("pool.value = ", pool.value))
outdir.value <- "outdir.value.input"
print(paste0("outdir.value = ", outdir.value))
R1trunclen <- R1trunclen.value.input
print(paste0("R1trunclen = ", R1trunclen))
R2trunclen <- R2trunclen.value.input
print(paste0("R2trunclen = ", R2trunclen))
minlen <- minlen.value.input
print(paste0("minlen = ", minlen))


if (concat.value == "true"){

    print("dada2 denoised reads will be concatenated, not merged")

} else {

    print("dada2 denoised reads will be merged, not concatenated")

}

print("Organizing file locations and sample names...")

# specify subdirectories
trimpath <- paste0(wdPath, paste0("/trimmed/",dataset.value, run.value, rc.value)) # directory containing fwd and rev primer-trimmed reads
print('Trimmed input files are named  \n')
list.files(trimpath)[1:20] # Make sure the trimmed files are in there

filtpath <- paste0(wdPath, paste0("/filtered/",dataset.value, run.value)) # directory containing fwd and rev quality-filtered reads
print('Filtered files are named \n')
list.files(filtpath)[1:20] # Should be empty at beginning of script

# This creates a list of file paths for the primer-trimmed input files
# It assumes forward and reverse fastq filenames have format:
# SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fastqFs.trimmed <- sort(list.files(trimpath, pattern="_R1.fastq.gz", full.names = TRUE)) # get forward reads (R1)
fastqRs.trimmed <- sort(list.files(trimpath, pattern="_R2.fastq.gz", full.names = TRUE)) # get reverse reads (R2)
if(length(fastqFs.trimmed) != length(fastqRs.trimmed)) stop("Forward and reverse files do not match.")
print('fastqFs.trimmed are in:\n')
list(fastqFs.trimmed[1:5])
list(fastqFs.trimmed[1:5])

# Rewrite file names for filtered reads:
# It replaces 'trimmed' with 'filtered' in the full path (subdirectory name) without changing the base file name
fastqFs.filtered <- fastqFs.trimmed
fastqRs.filtered <- fastqRs.trimmed
# Replace the file path name
fastqFs.filtered <- str_replace_all(fastqFs.filtered,'trimmed/','filtered/')
fastqRs.filtered <- str_replace_all(fastqRs.filtered,'trimmed/','filtered/')
print('fastqFs.filtered are in:\n')
list(fastqFs.filtered[1:5])
list(fastqFs.filtered[1:5])

# Extract simplified sample names
sample.names <- str_remove_all(basename(fastqFs.filtered),'_R1.fastq.gz')
names(fastqFs.filtered) <- sample.names
names(fastqRs.filtered) <- sample.names
print('Sample names are:\n')
print(sample.names[1:20])

####### Quality filtering + read truncation #######

print("Beginning pre-dada2 sequence filtering...")

# We have set out trunclen and minlen values earlier in the script

filterOutput <- filterAndTrim(fwd=fastqFs.trimmed, filt=fastqFs.filtered,
              rev=fastqRs.trimmed, filt.rev=fastqRs.filtered,truncLen=c(R1trunclen,R2trunclen),
              compress=TRUE, verbose=TRUE, multithread=TRUE, minLen = minlen, truncQ = 2)

saveRDS(filterOutput ,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "filterOutput", dataset.value,  run.value, ".RDS")))

print("Carrying out error rate determination...")

print("Read 1...")

####### Learn the Error Rates and denoise #######
errF <- learnErrors(fastqFs.filtered, multithread=TRUE)
saveRDS(errF,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "errF", dataset.value,  run.value, ".RDS"))) # save intermediate file

print(" Read 2...")

errR <- learnErrors(fastqRs.filtered, multithread=TRUE)
saveRDS(errR,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "errR", dataset.value,  run.value, ".RDS"))) # save intermediate file

print("Carrying out sequence dereplication, denoising and merging on a per sample basis...")

if (concat.value == "true"){

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(fastqFs.filtered[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(fastqRs.filtered[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR,
    justConcatenate = TRUE, verbose=TRUE)
    mergers[[sam]] <- merger
}
# Save merger file output
saveRDS(mergers,file=paste0(wdPath, paste("/dada2output/", outdir.value,"/", "mergers", dataset.value,  run.value, ".RDS",sep="")))

} else {

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(fastqFs.filtered[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(fastqRs.filtered[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR,
    minOverlap = 12, justConcatenate = FALSE, verbose=TRUE)
    mergers[[sam]] <- merger
}
# Save merger file output
saveRDS(mergers,file=paste0(wdPath, paste("/dada2output/", outdir.value,"/", "mergers", dataset.value,  run.value, ".RDS",sep="")))

}

print("Constructing sequence table...")

#######  Construct sequence table #######
seqtab <- makeSequenceTable(mergers)

print("Performing de novo chimera removal...")

#######  Remove Chimeras #######

# For my sequence data- the file is too big and this makes it computationally impossible to perform chimera detection... We will split it.

# Assume 'seqtab' is your sequence table
n_splits <- 50  # Number of splits
seqs <- colnames(seqtab)  # Extract sequence column names
seq_chunks <- split(seqs, cut(seq_along(seqs), n_splits, labels = FALSE))  # Split into 10 chunks

# Function to process each subset
process_chunk <- function(seq_subset) {
  if (length(seq_subset) == 0) return(NULL)  # Handle empty case
  sub_table <- mergers[, seq_subset, drop = FALSE]  # Extract chunk
  sub_table <- removeBimeraDenovo(sub_table, method = "consensus", multithread = TRUE, verbose = TRUE)  # Bimera removal
  return(sub_table)
}

# Apply function to all chunks
filtered_chunks <- lapply(seq_chunks, process_chunk)

# Combine the filtered tables
seqtab.nochim <- do.call(cbind, filtered_chunks)

####### Export nochim sequence table for subsequent merging
saveRDS(seqtab.nochim,file=paste0(wdPath, paste("/dada2output/", outdir.value,"/", "ASVtab.nochim", dataset.value,  run.value, ".RDS",sep="")))

print("Creating sequence tracking table...")

####### Report the number of reads that made it through each step in the pipeline: #######
getN <- function(x) {
  sum(getUniques(x))
}

track <- cbind(
  filterOutput,
  #sapply(dadaFs, getN), # not preserved in this version
  #sapply(dadaRs, getN), # not preserved in this version
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

colnames(track) <- c(
  "input",
  "filtered",
  #"denoisedF",
  #"denoisedR",
  "merged",
  "nonchim")
rownames(track) <- sample.names
head(track)

####### write outputs #######

print("Writing final dada2 asv table and fasta file outputs...")

# Give sequence headers more manageable names (ASV_1, ASV_2...)
# Sike we're not doing that
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste(">",colnames(seqtab.nochim), sep = "") # No sensible name just ASV codes
#asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# Usually we would carry out this code to give such manageable names, but since we are combining seqtab tables together, we do NOT want to do this until$
# for (i in 1:dim(seqtab.nochim)[2]) {
# asv_headers[i] <- paste(">ASV", i, sep="_")
# }

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "ASVs",dataset.value,  run.value,".fasta")))

#count table
asv_tab <- t(seqtab.nochim)
#row.names(asv_tab) <- sub(">", "", asv_headers)
head(asv_tab)

#add placeholder column, then row.names=FALSE
Taxa <- rownames(asv_tab)
head(Taxa)
data <- cbind(Taxa,asv_tab)
head(data)
write.table(
  data,
  file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "ASVtable", dataset.value,  run.value, ".csv")),
              sep=",", quote=F, row.names = FALSE)

write.table(
  data,
  file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "ASVtable", dataset.value,  run.value, ".tsv")),
              sep="\t", quote=F, row.names = FALSE)

# Tracking table
write.table(
track,
                file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "tracking_table", dataset.value,  run.value, ".csv")),
                            sep=",", quote=F, col.names=NA)

print("dada2.R complete")
