# Load packages and set working directory.

library(tidyverse)
library(jsonlite)
library(ivs)


setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//yeeI")

# Define function to read in HMMER files.
read_nhmmer = function(file){
  df = read.table(file, skip = 2, sep = "", colClasses = c("character", rep("NULL", 3), rep("character", 11), rep("NULL", 8)), fill = TRUE, row.names = NULL)
  df = df[,1:12]
  df = replace(df, df=='', NA)
  colnames(df) = c("accessions", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "strand", "eval", "score", "bias")
  df = drop_na(df)
  df = df %>% mutate_at(c("hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "eval", "score", "bias"), as.numeric)
  
  return(df)
}


# Upload lineages acquired from NCBI taxdump and taxonkit based on each assembly's taxid.
lineages = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//ref_lineage.txt", col.names = "taxID", header = FALSE) %>% 
  separate(taxID, into = c("organism", "domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "merge") %>% 
  separate(organism, into = c("taxID", "organism"), sep = "\\s", extra = "merge") %>%
  drop_na()

# Upload NCBI RefSeq Assembly database accessions (beginning with "GCF_") and corresponding NCBI Nucleotide (nuccore) database accessions (beginning with "NC_" or "NZ_".
df_GCF_nc = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//GCF_nuccore_reps_clean.csv")

# Upload assembly metadata for NCBI assemblies and select relevant data columns.
lines = readLines("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//prfB//Data//Bioinformatics//assembly_data_report.jsonl")
lines = lapply(lines, fromJSON)
lines = lapply(lines, unlist)
df_GCF_tax = bind_rows(lines)

df_GCF_tax_select = df_GCF_tax %>%
  select(accession, assemblyInfo.assemblyLevel, organism.organismName, assemblyInfo.bioprojectLineage.bioprojects.title, organism.taxId, checkmInfo.checkmMarkerSet, checkmInfo.completeness, checkmInfo.contamination, assemblyStats.totalSequenceLength, assemblyStats.gcPercent, assemblyStats.numberOfComponentSequences)

# Join the NCBI accessions with their metadata.
df_GCF_nc_tax = left_join(df_GCF_nc, df_GCF_tax_select, by = join_by(assembly == accession)) %>%
  mutate(checkmInfo.completeness = as.numeric(checkmInfo.completeness)) %>%
  mutate(checkmInfo.contamination = as.numeric(checkmInfo.contamination)) %>%
  mutate(assemblyStats.totalSequenceLength = as.numeric(assemblyStats.totalSequenceLength)) %>%
  mutate(assemblyStats.gcPercent = as.numeric(assemblyStats.gcPercent)) %>%
  mutate(assemblyStats.numberOfComponentSequences = as.numeric(assemblyStats.numberOfComponentSequences))

# Upload HMMER outputs and filter to remove proteins with low scores.
df_yeeI = read_nhmmer("yeeI_HMMER_5e2.csv") %>%
  mutate(start = ifelse(strand == "-", alito, alifrom)) %>%
  mutate(stop = ifelse(strand == "-", alifrom, alito)) %>%
  distinct(accessions, .keep_all = TRUE)

anti_yeeI = filter(df_yeeI, df_yeeI$score < 100)

ggplot(data = df_yeeI, aes(x = score)) +
  geom_histogram(bins = 37)

ggplot(data = df_yeeI_filt, aes(x = score)) +
  geom_histogram()

df_yrbC = read_nhmmer("yrbC_HMMER_5e2.csv") %>%
  mutate(start = ifelse(strand == "-", alito, alifrom)) %>%
  mutate(stop = ifelse(strand == "-", alifrom, alito)) %>%
  distinct(accessions, .keep_all = TRUE)

# Distinguish whether a shared gene hit is more like yeeI or yrbC and label it as such.
joined = full_join(df_yeeI, df_yrbC, by = "accessions", suffix = c("_yeeI", "_yrbC"))

check = joined %>%
  mutate(overlap = map2_lgl(iv(start_yeeI, stop_yeeI), iv(start_yrbC, stop_yrbC), type = "any", iv_overlaps)) %>%
  mutate(best_score = ifelse(overlap == "TRUE", pmax(score_yeeI, score_yrbC), NA))
check$best_score_gene = names(check[c("score_yeeI", "score_yrbC")])[max.col(check[c("score_yeeI", "score_yrbC")], "first")]
check$best_score_gene = sub("score_", "", check$best_score_gene)

df_filt = filter(check, score_yeeI > 120 | score_yrbC > 120)
df_yeeI_final = df_filt %>%
  filter(best_score_gene == "yeeI" | is.na(score_yrbC) | (overlap == FALSE & !is.na(score_yeeI)))

df_yrbC_final = df_filt %>%
  filter(best_score_gene == "yrbC" | is.na(score_yeeI) | (overlap == FALSE & !is.na(score_yrbC)))

# Prep df for extracting hit seqs.

select_yeeI = df_yeeI_final %>%
  select(accessions, ends_with("yeeI")) %>%
  rename_with(~str_remove(., '_yeeI'))

select_yrbC = df_yrbC_final %>%
  select(accessions, ends_with("yrbC")) %>%
  rename_with(~str_remove(., '_yrbC'))

df_bind = rbind(select_yeeI, select_yrbC)
#write_csv(df_bind, "yeeI_yrbC_hits_for_seq.csv")

# Match gene presence data with assembly and taxonomy data.
df = data.frame(cbind(df_GCF_nc_tax, yeeI = (df_GCF_nc_tax$nuccore %in% df_yeeI_final$accessions), yrbC = (df_GCF_nc_tax$nuccore %in% df_yrbC_final$accessions)))
  



# Multiple contigs, scaffolds, or chromosomes can make up an assembly. Identify gene hits across all contigs/scaffolds/chromosomes in each assembly. 
df_presence = df %>%
  group_by(assembly) %>%
  summarise(across(yeeI:yrbC, ~sum(.))) %>%
  mutate_if(is.numeric, ~1 * (. != 0))


bools = data.frame(cbind(assembly = df_presence$assembly, yeeI = ifelse(df_presence$yeeI == 1,"TRUE","FALSE"), yrbC = ifelse(df_presence$yrbC == 1,"TRUE","FALSE")))

# Save the names of all Nucleotide accession numbers that were searched in each assembly.
names =  df %>%
  distinct(nuccore, .keep_all = TRUE) %>%
  group_by(assembly) %>%
  summarize(nuccore=paste(nuccore, collapse=";"))

# Prepare the final dataframe with Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.
# Remove archaeal genomes. Filter based on CheckM completeness and contamination. Remove any duplicate genomes.
df_distinct = inner_join(names, bools, by = "assembly") %>%
  left_join(df[,1:11], by = "assembly", multiple = "any") %>%
  select(-X, -nuccore.y) %>%
  rename(nuccore = nuccore.x) %>% 
  inner_join(lineages, by = join_by("organism.taxId" == "taxID"), multiple = "all") %>%
  filter(!grepl("Archaea", domain)) %>%
  filter(checkmInfo.completeness > 80) %>%
  filter(checkmInfo.contamination < 10) %>%
  distinct(assembly, .keep_all = TRUE)

# Rename taxonomy to be consistent with Coleman et al 2021.
df_distinct$phylum[df_distinct$phylum == "delta/epsilon subdivisions"] = "Pseudomonadota"
df_distinct$phylum[df_distinct$phylum == "Pseudomonadota"] = "Proteobacteria"
df_distinct$phylum[df_distinct$phylum == "Terrabacteria group"] = df_distinct$class[df_distinct$phylum == "Terrabacteria group"]
df_distinct$phylum[df_distinct$phylum == "Bacillota"] = "Firmicutes"
df_distinct$phylum[df_distinct$phylum == "Actinomycetota"] = "Actinobacteriota"
df_distinct$phylum[df_distinct$phylum == "Abditibacteriota"] = "Armatimonadota"
df_distinct$phylum[df_distinct$phylum == "Aquificota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Campylobacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Deferribacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Thermodesulfobacteriota"] = "Desulfuromonadota + Desulfobacterota"
df_distinct$phylum[df_distinct$phylum == "Calditrichota"] = "FCB group"
df_distinct$phylum[df_distinct$phylum == "Thermomicrobiota"] = "Chloroflexota"
df_distinct$phylum[df_distinct$phylum == "Proteobacteria"] = df_distinct$class[df_distinct$phylum == "Proteobacteria"]

# Write Table S1.
write.csv(df_distinct, file = "yeeI_yrbC_df.csv")
