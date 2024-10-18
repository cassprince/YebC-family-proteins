library(ggprism)
library(tidyverse)
library(ggcorrplot)


df = read.csv("C:/Users/cassp/Box Sync/Feaga Lab/Cassidy Prince/yeeI/yebC_protein_identities_full.csv")

df2 = df %>%
  pivot_longer(cols = YebC_Ec:YeeI_Bs, names_to = "prot2", values_to = "identity") %>%
  separate(col = X, into = c("prot1", "org1")) %>%
  separate(col = prot2, into = c("prot2", "org2"))

df2$sub1 = paste0(df2$prot1, "[", df2$org1, "]")
df2$sub2 = paste0(df2$prot2, "[", df2$org2, "]")

#df2$sub1 <- factor(df2$sub1, levels=c("YebC[Ec]", "PmpR[Pa]", "YebC[Bs]", "YebC[Bb]", "YrbC[Bs]", "YebC[Ld]", "YeeN[Ec]", "YeeI[Bs]"))

df2$sub1 = with(df2, factor(sub1, levels = rev(unique(df2$sub1))))

ggplot(df2, aes(x = sub1, y = sub2, fill = identity)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="firebrick3", na.value = "white") +
  geom_text(aes(label = round(identity, digits = 2)), color = "black", size = 4) +
  theme_prism()+
  theme(text = element_text(size = 16))+
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
+
  scale_x_discrete(labels = parse(text = df2$sub2)) +
  scale_y_discrete(labels = parse(text = df2$sub2)) +
  labs(x = "", y = "")





ggsave("C:/Users/cassp/Box Sync/Feaga Lab/Cassidy Prince/yeeI/yebC_identities_matrix.png", p, width = 6, height = 4, units = "in")
