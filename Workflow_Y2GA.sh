#------------------------------------------------------------------------------------------
#---------------------------------- Required programs -------------------------------------
#------------------------------------------------------------------------------------------

# These programs are accessed throught the UGA GACRC resource high performance cluster (HPC)

#------------------------------------- HPC programs ---------------------------------------

# DSGSGSD

#-------------------------------------- R packaged ----------------------------------------

# sdgadg

#------------------------------------------------------------------------------------------
#------------------------------- HPC Directory Structure -------------------------------------
#------------------------------------------------------------------------------------------

# Almost all directories have sub-directories created as part of the pipeline based on user defined inputs
# I did this to account for the fact that sometimes I process multiple amplicon datasets, and it's nice to be able to run them all together in one master directory at the same time

# Directory: scripts
mkdir ./scripts # Holds all of the scripts written for this workflow

# Directory: trimmed
mkdir -p ./trimmed # Holds cutadapt trimmed sequences

# Directory: filtered
mkdir -p ./filtered # Holds quality filtered sequence outputs from the dada2 pipeline 'FilterandTrim' command

# q2files #Holds trained taxonomy classifier and outputs from 'asvfilter.sh' and 'taxonomy_qiime2.sh'
mkdir -p ./q2files

#
mkdir -p dada2output
# Sub-directories are made as part of...

mkdir -p readcounts

mkdir -p targetcuration

mkdir -p fastqcoutput

mkdir -p eukaryome # Holds various eukaryome databases downloaded through 'eukaryomedl.sh'

mkdir -p finaloutputs # Holds all final data relevant to downstream analysis

#------------------------------------------------------------------------------------------
#------------------------- Cutadapt trimming of primer sequences --------------------------
#------------------------------------------------------------------------------------------

# Sequencing was performed for all runs on an Illimina Nextseq 2000 (PE300)
# Amplicon fragments were frequently shorter than 300bp, resulting in through-reading by the sequencer
# To fix this we implement cutadapt twice:
#      - Once for removing primers in expected orientation from FWD and REV reads [scripts/cutadapt.sh]
#      - Again, this time removing the REVERSE COMPLEMENT of the OPPOSITE primer expected in the FWD and REV reads [scripts/cutadaptrc.sh]

# Determine number of samples for run 1
cd ./rawrun1/ # Move into the directory containing the raw sequences
ls -1 | wc -l
cd ../ # Move back into the working directoryinput

# Run cutadapt [run 1]
nSamples=1010
sbatch --export=fwdprimer=ACCCGCTGAACTTAAGC,revprimer=GTCGTTTAAAGCCATTACGTC,dataset=amflsu,run=run1 --array=1-$nSamples%10 ./scripts/cutadapt.sh

# Count the number of reads maintained after trimming
# cutadapt and cutadaptrc results should match- we did not discard reads in the second trimming
sbatch --export=dataset=amflsu,run=run1,RC= scripts/countreadsillumina.sh

#------------------------------------------------------------------------------------------
#------------------------------------ Quality Check ---------------------------------------
#------------------------------------------------------------------------------------------

# Let's see what the quality profiles of our sequences look like
# This script generates a DATASET-WIDE set of statistics, not individual sample statistics, it's for an overview
# We will use this to define if / where we perform sequence trimming in subsequent steps

# Carry out fastqc for all datasets
sbatch --export=dataset=amflsu,run=run1,RC= scripts/fastqc.sh

#------------------------------------------------------------------------------------------
#---------------------------- Error correction and denoising ------------------------------
#------------------------------------------------------------------------------------------

# We now generate amplicon sequence variants (ASVs) using the dada2 pipeline
# This covers multiple steps
#     - Quality control (trim sequences to defined length / by quality cut-off, remove Ns, filter seqs with Ns or below maxEE / maxN values.
#     - Dereplicate sequences (eases computational loads)
#     - Denoise sequences (error correction)
#     - Merge / Concatenate paired-end reads
#     - Remove chimeric sequences
#     - Construct sequence table
# We do this independly per sequence run because they may require different error models (based on per-run variability) for the best success

# dada2 denoising [run1]
sbatch --export=R1cutoff=220,R2cutoff=220,dataset=amflsurun1,RC=,minlen=220,directoryinput=amflsu,outdir=amflsu,concatenate=true ./scripts/dada2denoise_bigdata.sh

# Now let's use the mergeasvtabs file to do some otu table cleanup
sbatch --export=directoryinput=amflsu,dataset1=amflsurun1,dataset2=,dataset3=,dataset4=,runcount=1 scripts/mergeasvtabs.sh

#------------------------------------------------------------------------------------------
#---------------------------- Post-denoising quality control ------------------------------
#------------------------------------------------------------------------------------------

#--------------------------------- Frequency filtering ------------------------------------

# Lovely- we have sequences, but... are they sensible?
# We perform some light filtering to remove sequences that are unlikely to biologically correct. Especially important because we performed our dada2 processing in pooled mode (spurious seqs may show up as singletons in +1 samples)
# We filter by:
#     - otu abundance (>10 seqs in at least one sample)
#     - otu frequence (>0 seqs in at least 13 samples)
#     - sample filtering (samples > 1 sequence)
#           - In other cases I would set this as a realistic number, but we are going to merge multiple low seq depth samples for this study

sbatch --export=directoryinput=amflsu,dataset=amflsu,outdir=amflsu,otulevel=asv,otutotalfreqcutoff=10,otusamplefreqcutoff=13,minsamplefreqcutoff=1 ./scripts/asvfilter.sh

#------------------------------------ BLAST Filtering -------------------------------------

# Following reccomendations for removing non-AMF / biologically sensible sequences we perform BLAST filtering of sequences
# REFERENCE: doi.org/10.1111/nph.17080 and other Delavaux LSU database papers

# Filtering in this case is performed against the EUKARYOME LSU database from which I have subset glomeromycota sequence entries
# https://eukaryome.org
# REFERENCE: doi.org/10.1093/database/baae043 + doi.org/10.3897/mycokeys.107.125549

sbatch --export=targetdatabaseinput=./eukaryome/lsu/General_EUK_LSU_v1.9.3_glomeromycota.fasta,targetdatabaseblast=./eukaryome/lsu/General_EUK_LSU_v1.9.3_glomeromycota,outdir=amflsu,database=amflsu,otulevel=asv,dataset=amflsu,blastcuration=amfblast scripts/extractblast.sh

#------------------------------------------------------------------------------------------
#---------------------------------- Taxonomy Assignment -----------------------------------
#------------------------------------------------------------------------------------------

# We will assign taxonomy using a NAIVE BAYESIAN CLASSIFIER trained from the EUKARYOME LSU database, again subset for glomeromycota only
# There doesn't appear to be a pre-trained classifier available, so I trained one myself (train = TRUE), when I re-run this I use my pre-trained classifier (available in this repository).
#      - NOTE: IF YOU USE MY PRE-TRAINED CLASSIFIER- it is only relevant to the version of QIIME (and associated dependencies) that I have trained it on, if you use a different version, RETRAIN.
# My workflow for taxonomy:
#      - Assign taxonomy using naive bayesian classifier
#      - Filter out sequences that were NOT assigned to the kingdom 'Fungi' or phylum 'Glomeromycota'
#      - Perform taxa re-labeling (automated)
#             - Attach a unique 'ASV label identifier' for smooth downstream life e.g., 'Funneliformis mosseae ASV 1', 'Unclassified Rhizophagus ASV 4'
#             - For Kingdom-Species taxonomic levels, explicitely note where an ASV is 'Unclassified' at that level
#             - Neither of these are necessary but they sure do help with other filtering and interpretation downstream
#      - Decorate backbone phylogeny with AMF sequences
#             - RAxML tree generation based on previous tree made from LSU sequences (from EUKARYOME)
#      - Manually re-annotate GENUS / SPECIES level determinations WHERE APPROPRIATE with rationale

#---------------------------- Taxonomy assignment and filtering -----------------------------------

# Run naive bayesian classifier (QIIME2)
sbatch --export=outdir=amflsu,directoryinput=amflsu,otulevel=asv,dataset=amflsu,dbinput=eukaryome/lsu/QIIME2_EUK_LSU_v1.9.3.fasta,taxainput=eukaryome/lsu/QIIME2_EUK_LSU_v1.9.3_cleaned.tsv,curation=amfblast,training=FALSE scripts/taxonomy_qiime2.sh

# Perform filtering based on taxonomic assignment and relabeling
sbatch --export=kingdomfilter=Fungi,phylumfilter=Glomeromycota,directoryinput=amflsu,otulevel=asv,dataset=amflsu,curation=amfblast,refdbID=eukaryome,phylumunclassifiedfilter=FALSE ./scripts/taxonomyfilter_qiime2.sh

#------------------------- Align sequences against backbone alignment ------------------------------

# We have previously created a backbone phylogeny based on EUKARYOME lsu sequences
# See: Workflow_BackbonePhylogeny.sh for details
# I have copied this output into our working directory

# Align sequences against those used to construct the backbone phylogeny
sbatch --export=otulevel=asv,dataset=amflsu,backbonefile=eukaryome_backbonephylogeny/eukaryome/lsu/General_EUK_LSU_v1.9.3_glomeromycota_backbonesubset_headerreformat_outgroupcat.fasta,curation=amfblast scripts/alignseqsfinal.sh

#--------------------- Decorate backbone phylogeny with experimental sequences ---------------------

sbatch --export=dataset=amflsu,otulevel=asv,curation=,backbonetree=BACKBONE scripts/RAxMLtree.sh

#------------------------- Pull everything out and do some manual editing --------------------------

# Hate to do this to you but now comes the 'subjective' part of phylogeny-based taxonomic determination that doesn't require / allow devolvement to a reproducible code-based response...
# One day we'll have hundreds of representative seqs for each unique species / genus combination from which we can define edges to determine where our amplicon sequences fall between, today is not that day
# This is still better than the hundreds of spurious sequences defined as some kind of glomus in NCBI so let me live

# Visualise tree: iTol https://itol.embl.de
# I try to operate by some clear rules on deciding...