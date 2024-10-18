# Load packages and set working directory.
library(tidyverse)
library(treeio)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(castor)
library(ggbreak)
library(ggprism)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//yeeI")

# Define function to format tables for the tree.
table_for_tree = function(phylum, gene, accessions){
  table = as.data.frame(prop.table(table(phylum, gene), margin = 1)*100)
  table_new = data.frame(cbind(gene = table$Freq[table$gene == "TRUE"]))
  rownames(table_new) = accessions
  return(table_new)
}

# Upload dataframe containing NCBI Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.(Table S1).
df_GCF_nc_compressed = read.csv("yebC_family_hits_df.csv") 

# Separate into one Nucleotide accession number per row.
df_GCF_nc = separate_rows(df_GCF_nc_compressed, nuccore, sep = ";")

# Make a smaller dataframe with only data from phyla with >10 genomes.
dfShort = df_GCF_nc_compressed %>%
  group_by(phylum) %>%
  filter(n() > 10)

n_vals = dfShort %>%
  group_by(phylum) %>%
  summarize(n = n())

# Randomly sample one assembly per phylum for the collapsed tree. This was performed once and the assembly accessions were saved to the file "random_genomes_for_tree_new.csv" for reuse. 
random = dfShort %>%
  group_by(phylum) %>%
  sample_n(1) %>%
  select(-X) %>%
  ungroup()

rownames(random) = random$assembly
#write.csv(random, "random_genomes_for_tree_new.csv")

random = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//random_genomes_for_tree_new.csv")
random = column_to_rownames(random, 'X')

####### ---- Collapsed tree ---- #######

# Upload 16S tree of all assemblies and filter out any genomes included in the tree that were not surveyed in our study.
tree = read.newick("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//16S_fasttree.tre")
tree = get_subtree_with_tips(tree, only_tips = df_GCF_nc$nuccore)$subtree

# Replace Nucleotide accession tip labels with corresponding GCF accessions.
tree_tip_GCF = left_join(data.frame(tree$tip.label), df_GCF_nc, by = join_by("tree.tip.label" == "nuccore"), multiple = "any")
tree$tip.label = tree_tip_GCF$assembly

# Filter tree for the random representative genomes and midpoint root.
subtree_bar = get_subtree_with_tips(tree, only_tips = rownames(random))$subtree
tree_mid = midpoint(subtree_bar)

# Format tables for gene presence to append to the tree.
table_yeeI = table_for_tree(dfShort$phylum, dfShort$yeeI, random$assembly)
table_yrbC = table_for_tree(dfShort$phylum, dfShort$yrbC, random$assembly)

# Plot the collapsed tree and heatmap.
p = ggtree(tree_mid, size = 0.8) %<+% random + 
  xlim(0, 11.75) + 
  geom_tiplab(aes(label=phylum), align = TRUE, size = 4) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") + 
  new_scale_fill()

p1 = gheatmap(p, table_yeeI, offset = 6.7, width=0.7, font.size=3.5, colnames = TRUE, color=NA, colnames_angle = 90, colnames_offset_y = -0.2) + 
  scale_fill_viridis_c(option="A", direction = -1, name="Percent\nwith gene")


# Plot the full 16S tree.
tree_mid_big = midpoint(tree)

df_yeeI = data.frame(df_GCF_nc_compressed["yeeI"], row.names = df_GCF_nc_compressed$assembly) 
df_yeeI$yeeI = as.character(df_yeeI$yeeI)
df_yrbC = data.frame(df_GCF_nc_compressed["yrbC"], row.names = df_GCF_nc_compressed$assembly)
df_yrbC$yrbC = as.character(df_yrbC$yrbC)
df_phylum = data.frame(df_GCF_nc_compressed$phylum, row.names = df_GCF_nc_compressed$assembly)

p = ggtree(tree_mid_big, layout="circular", size=0.2) 

p1 = gheatmap(p, df_yeeI, offset=-0.2, width=0.125, font.size=1, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("0" = "white", "1" = "sienna1", "2" = "red4"), name="yebC2 gene count", na.value = "white") + 
  theme(text=element_text(size=18)) +
  new_scale_fill()

p2 = gheatmap(p1, df_yrbC, offset=0.375, width=0.125, font.size=1, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("0" = "white", "1" = "#C6DBEF", "2" = "#4292C6", "3" = "#08306B"), name="yebC gene count", na.value = "white") + 
  theme(text=element_text(size=18)) +
  new_scale_fill()


p2 = gheatmap(p1, df_phylum, offset=0.3, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  theme(text=element_text(size=15)) +
  scale_fill_discrete(name = "Phylum/group", na.value = "white")

ggsave("yebC_family_16S_tree_10_16_24.png", p2, units = "in", width = 17, height = 13, dpi = 600)

#### --- MISC DATA --- ###

test = df_GCF_nc_compressed %>%
  select(phylum, yeeI, yrbC) %>% 
  mutate_if(is.numeric, ~1 * (. != 0)) 

sum(test$yeeI | test$yrbC)


summary = df_GCF_nc_compressed %>% 
  group_by(phylum) %>%
  summarize(n = n(), yeeI = sum(yeeI), yrbC = sum(yrbC), yeeI_prop = round((sum(yeeI)/n())*100, 3), yrbC_prop = round((sum(yrbC)/n())*100, 3))

write_csv(summary, "yeeI_yrbC_by_phylum.csv")
