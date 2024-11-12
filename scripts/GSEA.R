# Set the working directory
setwd("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results")

# run in ther console
#options(future.globals.maxSize = 8000 * 1024^2)

# Install this libraries if you don't have it
library(dplyr)
library(Seurat)
library(patchwork)
library(stringi)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(ComplexHeatmap)


install.packages("msigdbr")

library(msigdbr)

# Visualize which species are present in the package
msigdbr_species()

all_genes <- msigdbr(species = "Sus scrofa")

head(all_genes)



# Replace NA gene symbols with corresponding Ensembl gene IDs
all_genes$gene_symbol[is.na(all_genes$gene_symbol)] <- all_genes$ensembl_gene[is.na(all_genes$gene_symbol)]

fgsea_sets <- all_genes %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets


markers <- read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table11.csv")

markers$X <- gsub("^Sus-Sscrofa-", "", markers$X)

ranks <- markers$avg_log2FC#sign(markers$avg_log2FC) * (-log10(markers$p_val)) #val_adj ,markers$avg_log2FC) * (-log10(markers$p_val)
names(ranks) <- markers$X

head(ranks)

plot(ranks)

# checking min and max
#max(ranks)
#min(ranks)
# fix min and max ranking
#max_ranking <- max(ranks[is.finite(ranks)])
#min_ranking <- min(ranks[is.finite(ranks)])

#ranks <- replace(ranks, ranks > max_ranking, max_ranking * 10)
#ranks <- replace(ranks, ranks < min_ranking, min_ranking * 10)

#ranks <- sort(ranks, decreasing = TRUE) # sort genes by ranking

#plot(ranks)

#max(ranks)
#min(ranks)


# run GSEA

GSEAres <- fgseaMultilevel(pathways = fgsea_sets, stats = ranks)#nproc = 4, scoreType = 'std'

# check result

head(GSEAres)

# top 6 enriched pathways
head(GSEAres[order(pval)])

sum(GSEAres$pval < 0.01, na.rm = T)#sum(GSEAres[pval < 0.01])

sum(GSEAres$padj < 0.01, na.rm = T)#sum(GSEAres[, padj < 0.01])


topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(fgsea_sets[topPathways], stats = ranks, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()




# only plot the top 20 pathways
ggplot(GSEAres %>% filter(padj < 0.01) %>% head(n= 30), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj < 0.01)) +#
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


# Dot plot of top pathways (p-value vs NES)
ggplot(GSEAres %>% filter(padj < 0.01) %>% head(n= 30), aes(reorder(pathway, NES), -log10(padj))) +
  geom_point(aes(size = NES, color = -log10(padj))) +  # Size by NES, color by adjusted p-value
  scale_color_gradient(low = "blue", high = "red") +  # Gradient color from blue (not significant) to red (high significance)
  labs(x = "Pathways", y = "-log10(Adjusted P-value)", 
       title = "Dotplot of GSEA Results (colored by padj)", 
       size = "NES", color = "-log10(padj)") +  # Add labels
  coord_flip() +  # Flip axes to make it horizontal
  theme_minimal()  # Clean theme


plotEnrichment(fgsea_sets[["BLANCO_MELO_BETA_INTERFERON_TREATED_BRONCHIAL_EPITHELIAL_CELLS_UP"]],
               ranks) + labs()

library(enrichplot)


GSEAres[GSEAres$pathway == "BLANCO_MELO_HUMAN_PARAINFLUENZA_VIRUS_3_INFECTION_A594_CELLS_UP"]$leadingEdge

# Read all the tables in a list

fgresList <- list(
  condition1 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table1.csv"),
  condition2 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table2.csv"),
  condition3 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table3.csv"),
  condition4 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table4.csv"),
  condition5 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table5.csv"),
  condition6 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table6.csv"),
  condition7 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table7.csv"),
  condition8 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table8.csv"),
  condition9 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table9.csv"),
  condition10 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table10.csv"),
  condition11 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table11.csv"),
  condition12 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table12.csv"),
  condition13 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table13.csv"),
  condition14 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table14.csv"),
  condition15 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table15.csv"),
  condition16 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table16.csv"),
  condition17 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table17.csv"),
  condition18 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table18.csv"),
  condition19 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table19.csv"),
  condition20 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table20.csv"),
  condition21 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table21.csv"),
  condition22 = read.csv("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/table22.csv")
)

# This function change the name of the genes for Scrofa and virus removing the analysis labelled added by cellranger
# Then it creates the rank using the log2FC and find the pathways

gsea_analisys <- function(markers, fgsea_sets) {
  print("Deleting subfix in front of the gene name")
  markers$X <- gsub("^Sus-Sscrofa-", "", markers$X)
  markers$X <- gsub("^h1n1--------", "", markers$X)
  print("Creating rakns of gene")
  ranks <- markers$avg_log2FC
  names(ranks) <- markers$X
  print("Calculatin GSEA result...")
  GSEAres <- fgseaMultilevel(pathways = fgsea_sets, stats = ranks)
  print("Done! Saving dataframe in a list")
  return(GSEAres)
}
# Apply it to all obeject in the list
gseare_list <- lapply(fgresList,gsea_analisys, fgsea_sets)

# Filter pathways found selecting the top 30 up/down regulated
gseare_list_padj <- lapply(gseare_list, function(res) {
  res <- res[res$padj < 0.01, ]  # Keep only pathways with padj < 0.01
  res <- res[abs(res$NES) > 2.0, ]
  # Select top 20 upregulated and top 20 downregulated pathways
  top_upregulated <- head(res[order(res$NES, decreasing = TRUE), ], 30)
  #print(top_upregulated)
  top_downregulated <- head(res[order(res$NES, decreasing = FALSE), ], 30)
  #print(top_upregulated)
  # Combine top pathways
  top_pathways <- rbind(top_upregulated, top_downregulated)
  #print(top_pathways)
 #print(res)
  return(res[res$pathway %in% top_pathways$pathway, ])
  #print(res)
  #return(res)
  
})

print(gseare_list_padj)



# Bar plot of top pathways (p-value vs NES)
ggplot(gseare_list_padj$condition22, aes(reorder(pathway, NES), NES)) +
  #geom_point(aes(size = size, color = padj)) +  # Size by NES, color by adjusted p-value
  geom_col(aes(fill = padj))+
  scale_fill_gradient(low = "red", high = "blue") +  # Gradient color from blue (not significant) to red (high significance)
  labs(x = "Pathways", y = "Normalized NES", 
       title = "Dotplot of GSEA Results", 
       size = "Count", color = "padj") +  # Add labels
  coord_flip() +  # Flip axes to make it horizontal
  theme_minimal()+  # Clean theme
  #scale_y_continuous(breaks = seq(-2,2,2))+
  theme(axis.text.y = element_text(size=5))

##################################################
print(gseare_list_padj)


# This is part can be changed according the interest of the data to visualize

# If you are interest in some condition to compare change the condition in this part
 gseare_list_filter <- list(
   condition11 = gseare_list_padj$condition11,
   condition12 = gseare_list_padj$condition12,
   condition15 = gseare_list_padj$condition15,
   condition16 = gseare_list_padj$condition16,
   condition17 = gseare_list_padj$condition17,
   condition18 = gseare_list_padj$condition18
)



# Get a list of all pathways in the GSEA results
all_pathways <- unique(unlist(lapply(gseare_list_filter, function(res) res$pathway)))

# Create a matrix to store NES values
nes_matrix <- matrix(NA, nrow = length(all_pathways), ncol = length(gseare_list_filter))
rownames(nes_matrix) <- all_pathways
colnames(nes_matrix) <- names(gseare_list_filter)

# Fill the matrix with NES values for each condition
for (i in 1:length(gseare_list_filter)) {
  res <- gseare_list_filter[[i]]
  nes_matrix[res$pathway, i] <- res$NES  # Fill in NES for matching pathways
}

# Convert NES matrix to a data frame for better handling (optional)
nes_df <- as.data.frame(nes_matrix)

# If you selected all the condition uncomment this part
# Create condition annotations
# annotation <- data.frame(Condition = c("WT Mock vs CF Mock","WT Mock 33°C vs CF Mock 33°C",
#                                        "WT Mock 37°C vs CF Mock 37°C",
#                                        "WT h1n1/h1n1-R38A vs CF h1n1/h1n1-R38A",
#                                        "WT h1n1-R38A vs CF h1n1",
#                                        "WT h1n1 vs CF h1n1-R38A",
#                                        "WT h1n1-R38A  33°C vs CF h1n1 33°C",
#                                        "WT h1n1-R38A 37°C vs CF h1n1 37°C",
#                                        "WT h1n1  33°C vs CF h1n1-R38A 33°C",
#                                        "WT h1n1 37°C vs CF h1n1-R38A 37°C",
#                                        "CF h1n1 vs CF Mock",
#                                        "CF h1n1-R38A vs CF Mock",
#                                        "CF h1n1 33°C vs CF Mock 33°C",
#                                        "CF h1n1 37°C vs CF Mock 37°C",
#                                        "CF h1n1-R38A 33°C vs CF Mock 33°C",
#                                        "CF h1n1-R38A 37°C vs CF Mock 37°C",
#                                        "WT h1n1 vs WT Mock",
#                                        "WT h1n1-R38A vs WT Mock",
#                                        "WT h1n1 33°C vs WT Mock 33°C",
#                                        "WT h1n1 37°C vs WT Mock 37°C",
#                                        "WT h1n1-R38A 33°C vs WT Mock 33°C",
#                                        "WT h1n1-R38A 37°C vs WT Mock 37°C"))

# If you selected specific condition, add here which one
# Create condition annotations
annotation <- data.frame(Condition = c("CF WT vs CF Mock","CF NS1 vs CF Mock","CF WT 33 vs CF Mock 33", "CF WT 37vs CF Mock 37","CF NS1 33vs CF Mock33","CF NS1 37vs CF Mock37"))
rownames(annotation) <- colnames(nes_matrix)


# change NA values in 0
nes_matrix[is.na(nes_matrix)] <- 0

# Plot all the condition all toghter
# If there are too many pathways it will be unreadable. Be aware!
# Use the ggplot (barplot) under
# Create the heatmap using NES scores (Uncomment)
# Heatmap(nes_matrix,
#         name = "NES",                        # Heatmap legend title
#         show_row_names = TRUE,               # Show pathway names
#         show_column_names = TRUE,            # Show condition names
#         cluster_rows = TRUE,                 # Cluster pathways
#         cluster_columns = TRUE,              # Cluster conditions
#         column_title = "Conditions",         # Title for the columns
#         row_title = "Pathways",              # Title for the rows
#         col = c("blue","white","red"),
#         top_annotation = HeatmapAnnotation(df=annotation),
#         row_names_gp = gpar(fontsize = 6),
#         heatmap_legend_param = list(title = "NES Score"))  # Customize legend
############################################################

# Dot plot of top pathways (p-value vs NES)
ggplot(nes_df, aes(reorder(pathway, condition1), -log10(condition1))) +
  geom_point(aes(size = condition1, color = -log10(condition1))) +  # Size by NES, color by adjusted p-value
  scale_color_gradient(low = "blue", high = "red") +  # Gradient color from blue (not significant) to red (high significance)
  labs(x = "Pathways", y = "-log10(Adjusted NES)", 
       title = "Dotplot of GSEA Results (colored by NES)", 
       size = "NES", color = "-log10(padj)") +  # Add labels
  coord_flip() +  # Flip axes to make it horizontal
  theme_minimal()  # Clean theme
nes_df$pathway <- rownames(nes_df)

