#load packages
library("Biostrings")
library(tidyverse)
library(seqinr)

setwd("/Users/lgschaer/Desktop/180_SIP_16S/")

fastaFile <- readDNAStringSet("steve_incubation.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
head(df)

barcodes <- read.table("steve_map.txt", col.names = c("SampleID",	"BarcodeSequence",	"LinkerPrimerSequence",	"SampleType",	"DepthCm",	"Incubation",	"Description"))
head(barcodes)

new_fasta <- df %>%
  separate(seq_name, into = c("SampleID", "SampleID2"), sep = "_") %>%
  full_join(barcodes) %>%
  mutate(
    Header = paste0(SampleID, "_", SampleID2, " orig_bc=", BarcodeSequence, " new_bc=", BarcodeSequence, " bc_diffs=0")
  ) %>%
  select(Header, sequence)
head(new_fasta)

dim(df)
dim(barcodes)
dim(new_fasta)

names <- as.character(new_fasta$Header)
sequences <- as.list(new_fasta$sequence)

write.fasta(sequences, names, "new_fasta.fa", open = "w", nbchar = 253, as.string = FALSE)

