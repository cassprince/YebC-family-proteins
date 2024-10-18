library(tidyverse)
library(seqinr)
library(readxl)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//yeeI")
file = "yebC_AA_100424.fasta"
df = read.csv("yeeI_yrbC_df.csv") %>%
  separate_rows(nuccore, sep = ";") %>%
  as.data.frame()

# MUST BE A "VALUE"
accs = df[, "nuccore"]

fasta = read.fasta(file, as.string = TRUE)
annot = getAnnot(fasta)
output = list()

# This will take forever. That's ok.
for (n in accs){
  output = append(output, fasta[grep(n, annot)])
}

write.fasta(output, names = names(output), file.out = "yebC_AA_filt_10_16_24.fasta")