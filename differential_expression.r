library(tidyverse)
library(vroom)
library(curl)
#BiocManager::install("DESeq2")
#install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz", repo=NULL, type="source")
library(DESeq2)
sessionInfo() 

count_table = read.table(file = "count_matrix.tsv", header = TRUE, sep = "\t")
count_matrix = as.data.frame.matrix(count_table) 
head(count_matrix)
count_matrix$X <- NULL
#count_matrix$X.1 <- NULL
condition = factor(c(rep("WT", 2), rep("SSS", 2)))
coldata = data.frame(row.names=colnames(count_matrix), condition)
head(coldata)
count_matrix <- na.omit(count_matrix)
dds = DESeqDataSetFromMatrix(countData = round(count_matrix), colData=coldata, design = ~condition)
dds = DESeq(dds)
res <- results(dds)

mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.01, na.rm = TRUE)

#Let's only select P-Value significant <0.05
resq0.05 = subset(res, padj<0.5)
#Or for genes that are regulated in one direction only
FC_UP_1.5 = subset(res,log2FoldChange>1.5)
FC_DOWN_1.5 = subset(res,log2FoldChange< -1.5)
#As 0.05 is an arbitrary cutoff, sometimes it's best to simply rank from most significant to least
res_ranked = res[order(res["padj"]),]

#We can now save our data in .csv format for probing in future
write.csv(res[order(res["padj"]),], file="Yeast_dRNA_WT_vs_dRNA_SSS_results.csv")

library(biomaRt)
listMarts()
ensembl = useMart("ensembl")
datasets = listDatasets(ensembl)
datasets
ensembl = useDataset("scerevisiae_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
filters[1:5,]

attributes = listAttributes(ensembl)

# Let's take a look at the list of attributes below.
#head(attributes)
annots_lookup = getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id","external_gene_name"),
  values = results_deseq$Ensembl_ID,
  mart = ensembl)
gene_counts_annots = left_join(results_deseq, annots_lookup, by = c("Ensembl_ID" = "ensembl_gene_id"))
head(gene_counts_annots)

vsd <- varianceStabilizingTransformation(dds, blind = T)
plotPCA(vsd)

#-------------Volcano plot

with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Saccharomyces cerevisiae WT vs SSS", 
               col = "gray", xlim=c(-3,3),ylim=c(0,8)))

# Now colour the points with genes that are significant padjusted<0.05

with(subset(resq0.05), points(log2FoldChange, -log10(padj), pch=20, col="blue"))

# Now colour the points with genes that are regulated >1.5-fold
with(subset(FC_UP_1.5), points(log2FoldChange, -log10(padj), pch=20, col="green"))

# Now colour the points with genes that are regulated >1.5-fold
with(subset(FC_DOWN_1.5), points(log2FoldChange, -log10(padj), pch=20, col="red"))

#Lets add some features

legend("topleft", legend=c("All genes", 
                           "Significant genes", 
                           "upregulated >1.5 fold",
                           "downregulated >1.5 fold"),
       col=c("gray", "blue", "green", "red"), lty=1:2, cex=0.8)
abline(v=1, col="red",lty=2)
abline(v=-1, col="red",lty=2)


#Gene Ontology -------------------
#Gene Ontology analysis using the ClusterProfileR package. But first, what exactly is gene ontology?
#An ontology is a formal representation of a body of knowledge within a given domain

library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)
#BiocManager::install("org.Sc.sgd.db")
library(org.Sc.sgd.db)

gGO <- groupGO(gene = as.character(rownames(resq0.05)),# Input
               OrgDb         = org.Sc.sgd.db,    # Organism
               keyType       = "ENSEMBL",       # Defining our input (Ensembl Gene IDs)
               level         = 3,               # The GO level  at which to bin our genes
               ont           = "BP",            # The domain of GO to use (i.e. MF, CC, BP)
               readable      = TRUE)            # If true, gene IDs are mapped to gene symbols

# We can now inspect the top 10 lines of our output
head(gGO)
# You will notice every GO term has a description, a number of genes that fall within that ontology, a ratio of this compared with the total input, and a full breakdown of these genes.

# We can also show this graphically.
barplot(gGO, showCategory=20)

#----------- Gene Ontology Over-representation
eGO <- enrichGO(gene          = rownames(resq0.05), # Input (sig genes)
                universe      = rownames(res),      # Our total list of genes
                OrgDb         = org.Mm.eg.db, # Organism
                keyType       = "REFSEQ", #"ENSEMBL", # Defining our input (Ensembl Gene IDs)
                ont           = "CC", # The domain of GO to use (i.e. MF, CC, BP)
                pAdjustMethod = "BH", # The method of multiple test correction (i.e. Benjamini Hochberg)
                pvalueCutoff  = 0.01, # The p-value threshold
                readable      = TRUE) # If true, gene IDs are mapped to gene symbols
head(eGO)
dotplot(eGO) + ggtitle("Dotplot for GO Over-representation Analysis")
