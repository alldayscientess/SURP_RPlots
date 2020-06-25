# This is my RScript for plotting in R
# We will be using a subselection of the RNA-seq data Jeff shared with everyone last time
# Thanks Jeff!
# gene.description_sub.csv

#########################################################
# In this session we are going to go over the following:
#   1. Review of how to read in and analyze data
#   2. Volcano Plots
#   3. Scatter Plots
#   4. Boxplots and t.tests
#   5. Heatmaps
########################################################

#   1. Review of how to read in and analyze data
# Start every session with setting your working directory

setwd("~/Downloads/")

# Then make sure to add the libraries you'll be needing for the commands in this Script
library(ggplot2)
library(gplots)
library(reshape2)

# Next we want to read in the data 
gene.transcript=read.csv("~/Downloads/gene.description_sub.csv")

#Check how many dimensions the dataset has (rows and columns)
dim(gene.transcript)

# Check what the column names are of your data
colnames(gene.transcript)

# Also check what class of data you are using
class(gene.transcript)

# To look at your data set you can use the function View()
View(gene.transcript)

# Remeber you can sort your data in increasing or decreasing order
View(gene.transcript[order(gene.transcript$BC_DOX_1_fpkm, decreasing=TRUE),])
View(gene.transcript[order(gene.transcript$GeneName, decreasing=FALSE),])


#   2. Volcano Plots
# The first plot we are going to make is a volcano plot to look at the significant genes
# in our dataset (Thanks Jeff!!!)
# Well walk through some usefull arguments to include when plotting

plot(x = gene.transcript$log2FoldChange.L8_DOXvsL8., 
     y = -log10(gene.transcript$pvalue.L8_DOXvsL8.), 
     las = 1, 
     xlim = c(-4,4),
     ylim = c(0,200),
     col = ifelse(gene.transcript$padj.L8_DOXvsL8<=0.05, "red", "black"),
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)",
     main = "Volcano Plot of L8")

legend("topright", 
       legend = c("significant"), 
       pt.bg = c("red"), 
       bty = "n",
       pch = 21, 
       cex = 0.8)

text(x = gene.transcript$log2FoldChange.L8_DOXvsL8., 
     y = -log10(gene.transcript$pvalue.L8_DOXvsL8.), 
     labels = gene.transcript$GeneName, 
     cex = 0.5, pos = 2)

###############################################################
# Exercise 1 
# Create a volcano plot using a differnt cell line
# Make sure to change the x and y axis to fit your data
# Change your colors and the shapes of your points
###############################################################

#   3. Scatter Plots
# In this section we are going do two things
# 1. Subselect our data
# 2. Convert our data frame in to a matrix

newdata=gene.transcript[,c(2,32:55)]
new_names=c("GeneName",rep("BC",3),rep("BC_DOX",3),rep("DL",3),rep("DL_DOX",3), 
            rep("HT",3),rep("HT_ZN",3), rep("L8",3),rep("L8_DOX",3))
colnames(newdata)<-new_names
colnames(newdata)
View(head(newdata))

# We want to turn the first column into the rownames of the dataset so that 
# we can convert the data frame into a matrix for future graphs

# First we need to use some nifty functions from the tidyverse
newdata2 = remove_rownames(newdata) %>% column_to_rownames(var="GeneName")

# Then we can convert our data frame to a matrix 
newdata3=as.matrix(newdata2)

# Confirm that the data is now a matrix

?????????

# We can look at the effect of drug on specific genes

data_colors2=c(rep("red",3),rep("blue",3), rep("green",3),rep("black", 3))
names2=c("BC","DL", "HT", "L8")

plot(x = newdata3["DPEP1",c(1:3,7:9,13:15,19:21)], 
     y = newdata3["DPEP1",c(4:6,10:12,16:18,22:24)], 
     col = data_colors2, 
     pch = 1, # you can change your shapes to whatever you like
     cex = 1.5,
     las = 1, 
     xlab = "Expression of DPEP1", 
     ylab = "Expression of DPEP1 + DOX",
     main = "Effect of DOX on DPEP1 Expression")
legend("topleft", legend = names2, pt.bg = data_colors2[c(1,4,7,10)], 
       pch= 21, cex = 1, horiz=T)

plot(x = newdata3["GAPDH",c(1:3,7:9,13:15,19:21)], 
     y = newdata3["GAPDH",c(4:6,10:12,16:18,22:24)], 
     col = data_colors2, 
     pch = 1, # you can change your shapes to whatever you like
     cex = 1.5,
     las = 1, 
     xlab = "Expression of GAPDH", 
     ylab = "Expression of GAPDH + DOX",
     main = "Effect of DOX on GAPDH Expression")
legend("topleft", legend = names2, pt.bg = data_colors2[c(1,4,7,10)], 
       pch= 21, cex = 1, horiz=T)

###############################################################
# Exersise 2 
# Create your own scatter plot using whatever gene you find 
# significant or interesting
# Use whatever colors or shapes you want!
###############################################################


# We can also look at the two variables (genes) using a scatter plot

# data_colors=c(rep("red",3),rep("pink",3), rep("blue",3), rep("lightblue",3),
#               rep("darkgreen",3),rep("lightgreen",3),rep("black", 3), rep("grey",3))
# names=c("BC", "BC_DOX", "DL", "DL_DOX", "HT", "HT_DOX", "L8", "L8_DOX")
# 
# plot(x = newdata3["GAPDH",],
#      y = newdata3["BRCA1",], 
#      col = data_colors, 
#      pch = 19,
#      las =1,
#      xlim = c(0,1500),
#      ylim = c(0,14),
#      xlab = "Expression of GAPDH", 
#      ylab = "Expression of BRCA1",
#      main = "Expression of BRCA vs. GAPDH")
# 
# legend("bottomleft", legend = names, pt.bg = data_colors[c(1,4,7,10,13,16,19,22)], 
#        pch= 21, cex = 0.7)


# text(x=newdata3["GAPDH",], y=newdata3["BRCA1",], labels = c(rep("BC",3), rep("BC_DOX",3), 
#                                                             rep("DL",3), rep("DL_DOX",3), 
#                                                             rep("HT",3),rep("HT_DOX",3), 
#                                                             rep("L8",3), rep("L8_DOX",3)), 
#      cex = 0.5, pos = 3)

# nd_100=head(newdata,100)
# nd_100=remove_rownames(nd_100)
# nd2=nd_100 %>% column_to_rownames(var="GeneName")
# rownames(nd2)


#   4. Boxplots and t.tests

# Let's make our first boxplot for just one cell line Untreated/Treated
boxplot(newdata3["ACTB",1:3], 
        newdata3["ACTB",4:6],
        col = c("purple", "gold"),
        names = c("BC_UT", "BC_DOX"),
        las =1)

# Now let's make boxplots for all the data sets
boxplot(newdata3["ACTB",1:3], 
        newdata3["ACTB",4:6], 
        newdata3["ACTB",7:9], 
        newdata3["ACTB",10:12], 
        newdata3["ACTB",13:15], 
        newdata3["ACTB",16:18], 
        newdata3["ACTB",19:21],
        newdata3["ACTB",22:24],
        ylim = c(1200,2400),
        names=c("BC", "BC_DOX", "DL", "DL_DOX", "HT", "HT_DOX", "L8", "L8_DOX"),  
        las=2,
        col = c("purple", "gold"),
        ylab = "Expression (FPKM)",
        main = "Expression of ACTB")
legend("topright", legend = c("Untreated", "Treated"), fill = c("purple", "gold"))


# That difference between BC and BC_DOX looks really significant
# We can no for sure by calculating a p-value using the t.test() function
t.test(newdata3["ACTB",1:3], newdata3["ACTB",4:6])

#We can store the output of the t-test into variables
newdata3_tt=t.test(newdata3["ACTB",1:3], newdata3["ACTB",4:6])

#We can call up specific elements of the t-test by using the $ 
newdata3_tt$p.value

#A fun little function for rounding numbers is to use the round() function
round(newdata3_tt$p.value, digits=4)

###############################################################
# Exersise 3 
# run a t.test on L8 and L8_DOX
# Make boxplots for L8 and L8_DOX
# In the title of the graph include the p-value in parenthasis

###############################################################

#   5. Heatmaps


heatmap.2(newdata3[1:10,1:6], 
          Rowv = FALSE, 
          Colv = FALSE, 
          labRow=FALSE, 
          labCol = FALSE,
          margins=c(3,3), 
          scale = "col", 
          col = hmcols_pink8, 
          keysize = 1.0,
          breaks = seq(-2, 2, length.out = 257), 
          density.info="none", 
          trace="none", 
          main = "test")

###############################################################
# Exersise 4 
# Make a heatmap on all the cell lines

###############################################################