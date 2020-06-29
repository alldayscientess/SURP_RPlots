# This is my RScript for plotting in R
# We will be using a subselection of the RNA-seq data Jeff shared with everyone last time
# Thanks Jeff!
# gene.description_sub.csv

#########################################################
# In this session we are going to go over the following:
#   1. Review of how to read in and analyze data
#   2. Volcano Plots
#   3. Save Volcano Plot to a PDF
########################################################

#   1. Review of how to read in and analyze data
# Start every session with setting your working directory

setwd("/Users/txc065/Downloads/SURP_RPlots-master") # Or however works for your computer

# Then make sure to add the libraries you'll be needing for the commands in this Script
library(ggplot2)
library(gplots)

# Next we want to read in the data 
# You need to include a path that goes to your files
# Uncomment if you want to include the entire path
#gene.transcript=read.csv("/Users/txc065/Downloads/SURP_RPlots-master/gene.description_sub.csv")

# Or if you're already in the same directory as your files you can just do
gene.transcript=read.csv("gene.description_sub.csv")

#Check how many dimensions the dataset has (rows and columns)
dim(gene.transcript)

# Check what the column names are of your data
colnames(gene.transcript)

# Also check what class of data you are using
class(gene.transcript)

# To look at your data set you can use the function View()
# Uncomment if you want to View() your data 

#View(gene.transcript)


#   2. Volcano Plots
# The first plot we are going to make is a volcano plot to look at the significant genes
# in our dataset (Thanks Jeff!!!)
# Well walk through some usefull arguments to include when plotting

# Plot your volcano plot first 
# including labeling red the genes that exceed an adjusted p-value of 0.05
plot(x = gene.transcript$log2FoldChange.L8_DOXvsL8., 
     y = -log10(gene.transcript$pvalue.L8_DOXvsL8.), 
     las = 1, 
     xlim = c(-4,4),
     ylim = c(0,200),
     col = ifelse(gene.transcript$padj.L8_DOXvsL8<=0.05, "red", "black"),
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)",
     main = "Volcano Plot of L8")

# Then add your labels to your plot only for genes whos
# adjusted p-values that are less than or equal to 0.05 
# AND 
# log2 fold change greater than 1.5 (or 2 -- whatever you think is most improtant)
text(x = gene.transcript$log2FoldChange.L8_DOXvsL8., 
     y = -log10(gene.transcript$pvalue.L8_DOXvsL8.), 
     labels = ifelse(gene.transcript$padj.L8_DOXvsL8. <= 0.05 & abs(gene.transcript$log2FoldChange.L8_DOXvsL8.) >= 1.5, 
                     as.vector(gene.transcript$GeneName), NA),
     cex = 0.5, pos = 3)

# Then add your legend! 
# (It doesn't matter which order you do the text/legend)
legend("topright", 
       legend = c("significant"), 
       pt.bg = c("red"), 
       bty = "n",
       pch = 21, 
       cex = 0.8)

#   3. Save Volcano Plot to a PDF

# First you add the function pdf() with your chonen filename and .pdf
pdf("L8_DOXvsL8_VolcanoPlot.pdf")

plot(x = gene.transcript$log2FoldChange.L8_DOXvsL8., 
     y = -log10(gene.transcript$pvalue.L8_DOXvsL8.), 
     las = 1, 
     xlim = c(-4,4),
     ylim = c(0,200),
     col = ifelse(gene.transcript$padj.L8_DOXvsL8<=0.05, "red", "black"),
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)",
     main = "Volcano Plot of L8")
text(x = gene.transcript$log2FoldChange.L8_DOXvsL8., 
     y = -log10(gene.transcript$pvalue.L8_DOXvsL8.), 
     labels = ifelse(gene.transcript$padj.L8_DOXvsL8. <= 0.05 & abs(gene.transcript$log2FoldChange.L8_DOXvsL8.) >= 1.5, 
                     as.vector(gene.transcript$GeneName), NA),
     cex = 0.5, pos = 3)
legend("topright", 
       legend = c("significant"), 
       pt.bg = c("red"), 
       bty = "n",
       pch = 21, 
       cex = 0.8)

# Then you close the plot with the function dev.off()
dev.off()