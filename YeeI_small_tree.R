library(dplyr)
library(tidyverse)
library(treeio)
library(ggtree)
library(readxl)
library(phangorn)

# Read in accessions
accessions <- read_excel("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/YeeI/accessions.xlsx")

#Read in tree
tree = read.tree("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/YeeI/yeei_fasttree.tre")
tree_mid = midpoint(tree)

#Change tip label names to species names
tree_mid$tip.label <- sub("\\.1.*", ".1", tree_mid$tip.label)
tree_mid$tip.label <- sub("^'", "", tree_mid$tip.label)
tree_mid$tip.label <- accessions %>%
  # Match the current tip labels (accessions) to the 'accession' column in df
  filter(accession %in% tree_mid$tip.label) %>%
  # Ensure the order matches the tree tip labels
  arrange(match(accession, tree_mid$tip.label)) %>%
  # Extract the 'name' column as the new tip labels
  pull(name)

tree_mid$node.label <- as.numeric(tree_mid$node.label)*100

# Plot tree 
ggtree(tree_mid) + 
  geom_tiplab(size=8) +
  geom_nodelab(nudge_x=0.04, size = 5) +
  geom_treescale(x=0, y=7.5, offset=0.1, linesize=1) +
  xlim(0,0.9)

# Save tree
# ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/YeeI/yeei_tree2.png")
