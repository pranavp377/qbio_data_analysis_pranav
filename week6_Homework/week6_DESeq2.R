#Exercise 1.1
BiocManager::install("DESeq2")
library("DESeq2")

BiocManager::install("maftools")
library(maftools)

library(TCGAbiolinks)

BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method

sum_exp <- GDCprepare(query)

#1.2
colData(sum_exp)$age_at_index #colData has patient information as columns
is.na(colData(sum_exp)$age_at_index)
counts = assays(sum_exp)$"HTSeq - Counts"
patient_data = colData(sum_exp)
#skip step 3
patient_data = patient_data[!is.na(patient_data$age_at_index), ]
patient_data
patient_data$barcode
counts = counts[ ,patient_data$barcode]
length(counts)

patient_data$age_category = ifelse(patient_data$age_at_index < 50, "young", "old")
patient_data$age_category = factor(patient_data$age_category, levels = c("young", "old"))

#1.3
if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
counts

counts_row_sums = rowSums(counts)
counts_row_sums
length(counts) #29376438
low_counts_mask = counts_row_sums >= 10
low_counts_mask
sum(!low_counts_mask) #4810 FALSE values
counts_true = counts[low_counts_mask, ]
counts_true

#2.1
dds = DESeqDataSetFromMatrix(countData = counts_true,
                             colData = patient_data,
                             design = ~age_category)
dds_obj = DESeq(dds)
resultsNames(dds_obj) 
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

#2.2
results
str(results)
head(results)

my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))

order_indices = order(my_df$y)
# we expect c(5, 1, 3, 2, 4) because:
# 1 is index 5
# 2 is index 1
# 3 is index 3
# 4 is index 2
# 5 is index 4
order_indices
my_df = my_df[order_indices, ]
my_df

#2.3
row_order = order(results$padj)
row_order
results = results[row_order, ]
results
head(results, n = 20)
#Gene RF00100 is significantly more expressed in old patients vs young patients, with a log2FoldChange of -5.45101. According to GeneCards, this is an RNA gene that is classified under RFam. There is not too much information on the gene.

#2.4
log2FoldChange_threshold = 1
padj_threshold = 0.05
log2FoldChange_mask = results$log2FoldChange > log2FoldChange_threshold
above_threshold_log = results[log2FoldChange_mask, ]
nrow(results)

padj_mask = results$padj < padj_threshold
below_threshold_padj = results[padj_mask, ]
subset_results = results[(above_threshold_log && padj_mask), ]

#2.5
fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

plot(x = results$log2FoldChange,
     y = -log10(results$padj),
     xlab = "log2 Fold Change (Young/Old)", 
     ylab = "-log10 Adjusted p-value",
     pch = 20) # smaller solid circles

abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")



##ggplot sample
fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05


library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in young",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")


volcano_plot

#2.6
write.csv(x = results,
          file = "/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/week6_Homework/results.csv",
          row.names = FALSE)
