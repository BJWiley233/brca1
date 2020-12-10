## first 10 samples
df <- data.frame(condition=c(rep("case", 5), rep("control", 5)),
                 SNP=c(1,1,1,0,0,1,0,1,1,1))
x=matrix(table(df), ncol=2)

chisq.test(x, simulate.p.value = T)
t.test(c(14,24,9,6,19), c(8,11,37,15,34))

class(table(df))
fisher.test(x)

## rest of the samples including first 10
## total 204
setwd("/home/coyote/JHU_Fall_2020/Genome_Analysis/Genome_project")
brca1.exons.vars <- read.csv("SRRs to test - All 204.csv")
pheno <- read.csv("SRRs to test - Pheno.csv")
SNP <- read.csv("SRRs to test - 41223094T C.csv", header = F)

library("dplyr")
new.df <- left_join(brca1.exons.vars, pheno, by = c("SRR.ID" = "Run"))

df <- new.df[, c("SRR.ID", "Variants.in.BRCA1.Exons", "CaseControl")]
snp.variant <- gsub("_.*", "", SNP$V1)

sum(df$SRR.ID %in% snp.variant)

df$SNP.T.C <- df$SRR.ID %in% snp.variant

## t-test for all variants in BRCA1
t.test(df$Variants.in.BRCA1.Exons[df$CaseControl=="Case"],
       df$Variants.in.BRCA1.Exons[df$CaseControl=="Control"])

## chi-square for true/false having SNP or not
table(df[, c("CaseControl", "SNP.T.C")])
chisq.test(table(df[, c("CaseControl", "SNP.T.C")]))
df$CaseControl
df$SNP.T.C

