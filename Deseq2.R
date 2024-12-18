library(BiocManager)

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(GEOquery)
library(tidyverse)

counts_data <- read.csv("Counts.csv", row.names = 1)

# read the sample info

sample_info <- read.csv("sample_info.csv",row.names = 1)

colnames(sample_info)

# making shore the row_name in sample_info matchs to coloum names in countsdata

all(colnames(counts_data) %in% rownames(sample_info))

# are they are in same order?

all(colnames(counts_data) == rownames(sample_info))

#set facter level

sample_info$Type <- factor(sample_info$Type)

levels(sample_info$Type)

sample_info$Type <- factor(sample_info$Type, levels = c("Effector", "Exhausted"))



dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = sample_info,
                              design = ~ Type) 


# Set the rederence for the type factor

dds$Type <- factor(dds$Type, levels = c("Exhausted","Effector"))

sample_info$Type

# keeping rows that have at least 10 read counts

keep <- rowSums(counts(dds))>= 10

dds <- dds[keep,]
dim(dds)


# step3 run DESeq2

dds <- DESeq(dds)

deseq_results <- results(dds)

# changing deseq object R object data frame

deseq_results <- as.data.frame(deseq_results)
class(deseq_results)
head(deseq_results)
names(deseq_results)

deseq_results$GeneName <- row.names(deseq_results)

names(deseq_results)
head(deseq_results)

deseq_results <- subset(deseq_results, select=c("GeneName","padj","pvalue","stat","lfcSE","log2FoldChange","baseMean"))
names(deseq_results)



write.table(deseq_results, file = "deseq_result_all.tsv", row.names = F , sep = "\t")


# extract the most differetially expressed genes 
#select genes with a significant change in gene expression(adjusted p value below 0.05)
# and log 2 fold chang 1< and >1

# step1 filter based on p adjusted value

filtered <- deseq_results %>% filter(deseq_results$padj < 0.05) 

# filtered based on fold change , head we will use  threshold of 1

filtered <- filtered %>% filter(abs(filtered$log2FoldChange) >= 1 )

dim(deseq_results)
dim(filtered)


# save deseq results ,save the originaldata and filterd one

write.csv(deseq_results ,"deseq_results_all.csv")
write.csv(filtered , "deseq_results_filtered.csv")

# order the results table by increasing p value

deseq_results_order <- deseq_results[order(deseq_results$padj),]

head(deseq_results_order)
write.csv(deseq_results_order , "deseq_filtered_ordered.csv", row.names = F, sep = ",")

deseq_results_order <- deseq_results[order(deseq_results$log2FoldChange),]
head(deseq_results_order)

write.csv(deseq_results_order , "deseq_filtered_ordered_log2FoldChange.csv", row.names = F, sep = ",")
# normalized counts



normalized_counts <- counts(dds ,normalize = TRUE)
head(normalized_counts)

deseq_results <- read.csv("deseq_results_all.csv")
filtered  <- read.csv( "deseq_results_filtered.csv")
write.csv(normalized_counts ,"normalized_counts.csv")



