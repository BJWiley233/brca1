library(biomaRt)

human.mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl",
                      host = "useast.ensembl.org")

dat1 <- getBM(attributes = c("hgnc_symbol", "percentage_gene_gc_content", "scanprosite"),
              filters = c("chromosome_name", "band_start", "band_end"),
              values = list(17, "q21.31", "q21.31"),
              mart = human.mart,
              uniqueRows = T)

library(dplyr)

unique.genes <- dat1 %>% 
  select(hgnc_symbol, percentage_gene_gc_content) %>%
  distinct(hgnc_symbol, percentage_gene_gc_content) %>%
  arrange(desc(percentage_gene_gc_content))

brca1.gc <- unique.genes[unique.genes$hgnc_symbol=="BRCA1", "percentage_gene_gc_content"]

# BRCA1 percentile
sum(unique.genes$percentage_gene_gc_content <= brca1.gc)/nrow(unique.genes)

# Top 10 GC content
topGenes <- dat1 %>% 
  select(hgnc_symbol, percentage_gene_gc_content) %>%
  distinct(hgnc_symbol, percentage_gene_gc_content) %>%
  arrange(desc(percentage_gene_gc_content)) %>%
  top_n(10)

# Prosite Profiles
dat1 %>% 
  select(hgnc_symbol, scanprosite) %>%
  filter(scanprosite != "") %>%
  distinct()


#simple plot

library(ggplot2)


ggplot(unique.genes, aes(hgnc_symbol, percentage_gene_gc_content,
                         fill=factor(ifelse(hgnc_symbol=="BRCA1", "highlighted", "normal")))) +
  geom_histogram(stat="identity") +
  ggtitle("GC Content of Genes of Band 17q21.31") +
  labs(x="HGNC Symbol", y="GC Content") +
  theme(axis.text.x=element_text(angle=90),
        legend.position="bottom",
        plot.title=element_text(hjust=0.5, size=22))
