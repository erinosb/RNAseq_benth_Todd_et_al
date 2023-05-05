#######################################

#__Date:__ December
#__Author:__ Marc Nishimura / Erin Osborne Nishimura
#__Script:__ 240426_HopBA1_Deseq2.R
#__Project:__ To analyze RNA-seq in a project that compares Nicotiana benthiana treated with either HopBA1 or empty vector. Both leaf and petiole samples are included

#__Requires:__ 
# 
# + R (4.2.2)
# + DESeq2 (1.38.1)   
# + corrplot (0.92)
# + RColorBrewer (1.1-3)
# + pheatmap (1.0.12)
# + apeglm (1.20.0)
# + ashr () #this one is new

######################################


######### FOR FIRST TIME USE ONLY ##############
######### After use, comment this section ##############

### If you don't have bioconductor, install it. This version works with R version 4.1. Use "3.11" for R version 4.0
#if (!require("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install(version = "3.16")


### Install required packages/libraries:

### Install DESeq2:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")

#BiocManager::install("DESeq2")


### Install 'apeglm'
#BiocManager::install("apeglm")

### install corrplot:
#install.packages("corrplot")

###Install pretty heatmap - pheatmap: https://cran.r-project.org/web/packages/pheatmap/index.html
#install.packages("pheatmap")

###Install RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
#install.packages("RColorBrewer")

### SPLITCONDITIONS: Install ashr
install.packages("ashr")


######### DONE WITH: FIRST TIME USE ONLY SECTION ########################
######### ONCE YOU HAVE COMPLETED THAT SECTION ONCE, COMMENT IT OUT #####


###########  LOAD PACKAGE  #####################

# Do this every time.
# Load packages:
library(DESeq2)
library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(apeglm)

### SPLITCONDITIONS: Load ashr
library("ashr")
#################################################




###########  READ IN THE DATA  #####################

# import the counts data
getwd()

#Set this to your working directory:
# You may need to set this to your own working directory to your scripts directory:
#setwd("/set/to/your/scripts/directory")
#setwd("~/Dropbox/COURSES/2022_DSCI512/07_testrun_DESEq/02_scripts")
#setwd("C:/Users/MarcNishimura/Documents/DEseq2_HopBA1_Nb_20230418/02_scripts")

getwd()
countsData <- read.table(file = "../01_input/20230421_counts.txt", header = FALSE, row.names = 1, skip = 2) # 

# :!: EXERCISE

# Explore the countsData object using head and dim and class

head(countsData)
dim(countsData)
class(countsData)

# Read in the metadata
metadata1 <- read.table(file = "../01_input/HopBA1_Nb_RNAseq_metadata_20230415.txt", header = FALSE) # import the data
metadata1


# Organize the metadata file by adding somme column headers
colnames(metadata1) <- c("fasta1", "fasta2", "names1", "names2", "treatment", "tissue", "rep")
metadata1

# Exercise: Looking at the metadata1 object and countsData objects (try head, dim, str, class, colnames):


# Organize the countsData file.
# Notice that the countsData file doesn't have any column headers:
head(countsData)

# Let's give countsData some columns names. The first names will be... chr', 'start', etc...
# The last names will be names for each sample. We can pull those names from metadata1:
as.vector(metadata1$names2)



# Name countsData columns headers:
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata1$names2))

# :!: EXERCISE: Now look at the top of countsData again usign head():

head(countsData)



################### COUNT MATRIX INPUT ###################

# In this section we will prepare our input data for analysis.
# In the instruction for DESeq2, it states: "We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata."

# OK, our task will be to generate a table called "cts" out of the countsData table.
# Subset the countsData 
head(countsData)
dim(countsData)
#head(countsData[,6:23])
#dim(countsData[,6:23])
head(countsData[,6:17])
dim(countsData[,6:17])

# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:17])
head(cts)

# Yay, we made cts

# Next we need to make an information called coltable. We can make this out of the metadata table.

class(metadata1)
# Reorganize the metadata table so the names2 column are now row headers
metadata1
rownames(metadata1)<- metadata1$names2
metadata1$names2

# SPLITCONDITIONS: Add an additional column to allow us to split the condition groups later
metadata1 <- cbind(metadata1, combocond = as.factor(c("HopBA1_L", "HopBA1_L", "HopBA1_L", "HopBA1_P", "HopBA1_P", "HopBA1_P", "EV_L", "EV_L", "EV_L", "EV_P", "EV_P", "EV_P")))
metadata1$combocond <- factor(metadata1$combocond, levels = c("EV_L", "HopBA1_L", "EV_P", "HopBA1_P"))

# SPLITCONDITIONS: add in the combocond column to the coldata object
coldata <- metadata1[,c("treatment", "tissue", "combocond", "rep")]
coldata$treatment <- as.factor(coldata$treatment)
coldata$tissue <- as.factor(coldata$tissue)
coldata$rep <- as.factor(coldata$rep)

str(coldata)



# Yay! Now we have coldata! This is a new metadata object where we have just selected the type of information that is critical for deseq2 to use.

# One thing we need to explicitly check. The rownames of coldata need to exactly match the colnames of cts.
rownames(coldata)
colnames(cts)

#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))


# Next we will create an ddsHTSeq object out of cts and coldata:
# This will set a base design for your experiment:
# Load all the _counts.txt files and to attach them to the metadata.

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment + tissue)


################# PRE-FILTERING -- FILTER FOR PRESENT GENES: ################# 
# Not necessary, but helps keep things fast. 
# Exclude all samples that have less than 10 reads:
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Exercise: How many did we exclude?
dim(dds)

str(dds)

################### NOTE ON FACTOR LEVELS ###################
# Organize the categories into an order based on what makes sense:
dds$treatment <- factor(dds$treatment, levels = c("HopBA1","EV"))
dds$tissue <- factor(dds$tissue, levels = c("L", "P"))



# PERFORM DESEQ2 analysis:
####### This is where the magic happens! ########### 
# This will transform dds into a specialized object with many more fields filled in.
dds <- DESeq(dds)



# Exercise: Check the output (class, str, dim, plotDispEsts) by adding the object dds inside the parentheses:
#Write out your own code here:
class(dds)
str(dds)
dim(dds)
plotDispEsts(dds)


# Here is a demonstration of the size Factor scaling that was calculated (sizeFactor):
dds$sizeFactor

## NOTE:  09_EV_L_3 sample is an outlier here - 0.3156 



# Demo: Access the normalized counts using counts(x, normalized = TRUE)
# Demo: Access the raw count info using counts(x, normalized = FALSE)
# Both of these datasets are good supplemental tables for papers
head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))



############## DIFFERENTIAL EXPRESSION ANALYSIS #####################

# calculate the statistically significant differences between E.coli and B. subtilis diets
res_HopBA1_v_EV <- results(dds,
                            lfc = 0.5,
                            contrast=c("treatment", "EV", "HopBA1"))

resultsNames(dds)

summary(res_HopBA1_v_EV)


############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################

# An input requirement of the lfcShrink function is a coef term. This is pulled from the resultsNames of dds:

resultsNames(dds)

reslfc_HopBA1_v_EV <- lfcShrink(dds = dds, lfc = 0.5, coef="treatment_EV_vs_HopBA1", res = res_HopBA1_v_EV, svalue = FALSE)
help(lfcShrink)
summary(res_HopBA1_v_EV)
summary(reslfc_HopBA1_v_EV)
reslfc_HopBA1_v_EV



##################  Exploring and exporting results ##################  

##### KNOWN GENES:


##### MA PLOTS:

# Plot the the defaul MA-plot and the shrunken MA-plot:

par(mfrow=c(1,1))
plotMA(res_HopBA1_v_EV, main="HopBA1 v. EV \nunshrunken", ylim = c(-8,8), 
       ylab = "log fold change (ratio of normalized HopBA1 / EV)",
       xlab = "means of normalized counts")

plotMA(reslfc_HopBA1_v_EV, main="HopBA1 v. EV \nshrunken", ylim = c(-8,8), 
       ylab = "log fold change (ratio of normalized HopBA1 / EV)",
       xlab = "means of normalized counts")


# Save the plots

# Make a folder ../03_output
pdf("../03_output/MAplots.pdf", height = 6, width = 11)

par(mfrow=c(1,2)) # This changes how many plot panels you can generate
plotMA(res_HopBA1_v_EV, main="HopBA1 v. EV \nunshrunken", ylim = c(-8,8), 
       ylab = "log fold change (ratio of normalized HopBA1 / EV)",
       xlab = "means of normalized counts")

plotMA(reslfc_HopBA1_v_EV, main="HopBA1 v. EV \nshrunken", ylim = c(-8,8), 
       ylab = "log fold change (ratio of normalized HopBA1 / EV)",
       xlab = "means of normalized counts")


dev.off() 
dev.off() #Note: sometimes you need 2x dev.off() lines of code to really truly escape out of the .pdf printer

##### CORRELATION MATRICES:

#Take r-stabilized log transformations of all the normalized count data. This will help with the problem that the data is noisy and it will help with the problem that the data is spread across a wide range of values.
rld <- rlog(dds, blind=FALSE)  #Take the r-stabilized log transformed data:

# Calculate the distances between each sample
sampleDists <- dist(t(assay(rld))) # calculate distance matrices:


sampleDistMatrix <- as.matrix(sampleDists) #convert from data.frame -> matrix
rownames(sampleDistMatrix) <- colnames(rld) # Add some labels
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #Pick some pretty colors

# Draw the heatmap
par(mfrow=c(1,1))
p <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors) # Plot the heatmap

# Save the CORRELATION MATRIX as a .pdf file


pdf("../03_output/corr_matrix_plots.pdf", height = 6, width = 7)
        par(mfrow=c(1,1))
        p

dev.off() # Do you see "null device"
dev.off() #-> sometimes you need a second dev.off until you see "null device"


##### PCA PLOT ########################

# Draw a PCA Plot
par(mfrow=c(1,1))
plotPCA(rld, intgroup=c("treatment", "tissue"))

# Save a PCA Plot to a file
pdf("../03_output/PCAplot.pdf", width = 6, height = 6)
plotPCA(rld, intgroup=c("treatment", "tissue"))

dev.off() 
dev.off()

##### VOLCANO PLOTS:###################

resultsNames(dds)
# Volcano plots are nice ways of displaying the fold change against the p-value.
res_HopBA1_v_EV

significantLFC <- rbind(subset(res_HopBA1_v_EV, padj < 0.01 & log2FoldChange > 1.0  ), subset(res_HopBA1_v_EV, padj < 0.01 & log2FoldChange < -1.0 ))# Identify significantly changing genes

significant_points_to_plot <-res_HopBA1_v_EV[which(rownames(res_HopBA1_v_EV) %in% rownames(significantLFC)),] 

# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout <- subset(res_HopBA1_v_EV, padj < 10e-100)

#Draw plot:
par(mfrow=c(1,1)) # one plot only 

# Draw the plot
plot(res_HopBA1_v_EV$log2FoldChange, -log10(res_HopBA1_v_EV$padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-7,7), col = "#00000030")

# Add points
points(maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot$log2FoldChange, -log10(significant_points_to_plot$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")

# Add lines
abline(v=0, col = "blue")
abline(v=1.0, col = "blue", lty = "dashed")
abline(v=-1.0, col = "blue", lty = "dashed")
abline(h=-log10(0.01), col = "blue", lty = "dashed")



## Save the VOLCANO PLOT as a .pdf
pdf("../03_output/volcanoplot.pdf", width = 6, height = 6)

        par(mfrow=c(1,1)) # one plot only 

        # Draw the plot
        plot(res_HopBA1_v_EV$log2FoldChange, 
             -log10(res_HopBA1_v_EV$padj),
             main="Volcano plot", 
             xlab="Effect size: log2(fold-change)", 
             ylab="-log10(adjusted p-value)", 
             pch=20, 
             cex = 0.4, 
             ylim = c(0, 100), 
             xlim = c(-8,8), 
             col = "#00000030")

        # Add points
        points(maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)), 
               pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

        points(significant_points_to_plot$log2FoldChange, 
               -log10(significant_points_to_plot$padj), pch=20, cex = 0.4, 
               ylim = c(0, 100), col = "red")

        # Add lines
        abline(v=0, col = "blue")
        abline(v=0.01, col = "blue", lty = "dashed")
        abline(v=-0.01, col = "blue", lty = "dashed")
        #MTN what is this -log10 line?
        abline(h=-log10(0.1), col = "blue", lty = "dashed")


dev.off() 
dev.off()



############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################

# We will use the set of significantly changing genes from our variance shrunken analysis:
#Remember we calculated this above as:
#res_EcolVBsubt <- results(dds,
#                          lfc = 0.5,
#                          contrast=c("food", "Ecoli", "Bsubtilis"))
# resLFC_EcolVBsubt <- lfcShrink(dds, coef="food_Ecoli_vs_Bsubtilis", res = res_EcolVBsubt, type='apeglm')

res_HopBA1_v_EV
reslfc_HopBA1_v_EV

# Check the results table:
summary(reslfc_HopBA1_v_EV)
head(reslfc_HopBA1_v_EV)
str(res_HopBA1_v_EV)



# Select the significant subset of genes that are up-regulated in acetic acid
Up_in_A_1 <- subset(res_HopBA1_v_EV, padj < 0.01 & log2FoldChange > 1.0)
Up_in_A_1
dim(Up_in_A_1)
Up_in_A_1 <- Up_in_A_1[order(Up_in_A_1$log2FoldChange, decreasing = TRUE),] #order them
write.table(Up_in_A_1, file = "240412_MarcsList_1xlog2.txt", quote = FALSE)

# Select the significant subset of genes that are up-regulated in acetic acid
Up_in_A_2 <- subset(res_HopBA1_v_EV, padj < 0.01 & log2FoldChange > 2.0)
Up_in_A_2 <- Up_in_A_2[order(Up_in_A_2$log2FoldChange, decreasing = TRUE),] #order them
dim(Up_in_A_2)
write.table(Up_in_A_2, file = "240412_MarcsList_2xlog2.txt", quote = FALSE)


# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_A_1 <- subset(res_HopBA1_v_EV, padj < 0.01 & log2FoldChange < -1.0)
Down_in_A_1 <- Down_in_A_1[order(Down_in_A_1$log2FoldChange, decreasing = TRUE),]
dim(Down_in_A_1)


# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_A_2 <- subset(res_HopBA1_v_EV, padj < 0.01 & log2FoldChange < -2.0)
Down_in_A_2 <- Down_in_A_2[order(Down_in_A_2$log2FoldChange, decreasing = TRUE),]

head(Down_in_A_2)
dim(Down_in_A_2)

# Save these lists to output files:
write(rownames(Up_in_e.coli_2), file = "../03_output/Up2x_HopBA1.txt", sep = "\n")
write(rownames(Up_in_e.coli_1), file = "../03_output/Up1x_HopBA1.txt", sep = "\n")
write(rownames(Down_in_e.coli_1), file = "../03_output/Down1x_HopBA1.txt", sep = "\n")
write(rownames(Down_in_e.coli_2), file = "../03_output/Down2x_HopBA1.txt", sep = "\n")
write(rownames(resLFC_EcolVBsubt), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.

# Fun thing to do. See what KEGG pathways or GO Ontology terms are associated with your different lists of genes.
# Go to DAVID: https://david.ncifcrf.gov/
#      Click on "Start Analysis
#      Copy and paste your WBGENEID list into the "Upload" tab and "A: Paste a list" field
#      Under "Step 2: Select Identifier" select "WORMBASE_GENE_ID"
#      Under "Step 3:" select "Gene List"
#      Submit

######## SPLITCONDITIONS: Create Combinations of Conditions #####

# SPLITCONDITIONS - save dds object and re-specify the conditions

CCdds <- dds

# SPLITCONDITIONS - check the levels
levels(CCdds$treatment)
levels(CCdds$tissue)
levels(CCdds$combocond)

# SPLITCONDITIONS - specify a new type of condition
# SPLITCONDITIONS - this will create pairwise comparisons between the conditions specified in combocond
design(CCdds) <- formula(~ combocond)
CCdds <- DESeq(CCdds)


# SPLITCONDITIONS - Now we can perform contrasts of HopBA1 v. EV within leaves:
results_in_leaves <- results(CCdds, lfc = 0.5, contrast=c("combocond", "EV_L", "HopBA1_L"))
summary(results_in_leaves)

reslfc_in_leaves <- lfcShrink(dds = CCdds, lfc = 0.5, contrast=c("combocond", "EV_L", "HopBA1_L"), res = results_in_leaves, type="ashr")
summary(reslfc_in_leaves)

# SPLITCONDITIONS - and in petioles:
results_in_petioles <- results(CCdds, lfc = 0.5, contrast=c("combocond", "EV_P", "HopBA1_P"))
summary(results_in_petioles)

reslfc_in_petioles <- lfcShrink(dds = CCdds, lfc = 0.5, contrast=c("combocond", "EV_P", "HopBA1_P"), res = results_in_petioles, type="ashr")
summary(reslfc_in_petioles)

# SPLITCONDITIONS - MA-plots

pdf("../03_output/MAplots_for_leaves_and_petioles.pdf", height = 6, width = 11)

par(mfrow=c(1,2)) # This changes how many plot panels you can generate
plotMA(reslfc_in_leaves, main="HopBA1 v. EV shrunken \n within leaves", ylim = c(-8,8), 
       ylab = "log fold change (ratio of normalized HopBA1 / EV)",
       xlab = "means of normalized counts")

plotMA(reslfc_in_petioles, main="HopBA1 v. EV shrunken \n within petioles", ylim = c(-8,8), 
       ylab = "log fold change (ratio of normalized HopBA1 / EV)",
       xlab = "means of normalized counts")


dev.off() 
dev.off() #Note: sometimes you need 2x dev.off() lines of code to really truly escape out of the .pdf printer


######## SPLITCONDITIONS: Create a heatmap #####


# Let's loosen our restrictions on significance to all genes with any log fold change and adjusted p-values less than 0.1 (both are default)

# Get leaf-specific, HopBA1 v. EV differentially expressed genes:
results_in_leaves_2 <- results(CCdds, contrast=c("combocond", "EV_L", "HopBA1_L"), lfc =1)
summary(results_in_leaves_2)
# Get petiole-specific, HopBA1 v. EV differentially  expressed genes:
results_in_petioles_2 <- results(CCdds, contrast=c("combocond", "EV_P", "HopBA1_P"), lfc =1)
summary(results_in_petioles_2)
# Get differentially expressed genes between leaf v. petiole for empty vector treatments:
results_EV_2 <- results(CCdds, contrast=c("combocond", "EV_P", "EV_L"), lfc =1)

# Get differentially expressed genes between leaf v. petiole for HopBA1 treatments:
results_HopBA1_2 <- results(CCdds, contrast=c("combocond", "HopBA1_P", "HopBA1_L"), lfc =1)
summary(results_HopBA1_2)
#Subset each results table for just the differentially expressed genes:
sign_leaf <- subset(results_in_leaves_2, padj < 0.001)
sign_petiole <- subset(results_in_petioles_2, padj < 0.001)
sign_EV <- subset(results_EV_2, padj < 0.001)
sign_Hop <- subset(results_HopBA1_2, padj < 0.001)

dim(sign_leaf)
dim(sign_petiole)
dim(sign_EV)
dim(sign_Hop)

#Determine how many genes were captured and merge them:
changing_genes <- rbind(sign_leaf, sign_petiole)

dim(changing_genes)
length(unique(rownames(changing_genes)))
# Get all r-stabilized log values of all the data:
rdl_all <- assay(rlog(CCdds, blind=FALSE))

#Subset just the changing genes:
changing_lrt_rdl <- subset(rdl_all, rownames(rdl_all) %in% rownames(changing_genes))
dim(changing_lrt_rdl)
head(changing_lrt_rdl)

# Make sure it is in matrix form:
class(changing_lrt_rdl)
text1 = dim(changing_lrt_rdl)[1]

# Draw a heat map #1 - 
# Scale by row, listed below as scale = "row". This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:
plot_q <- pheatmap(changing_lrt_rdl, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=TRUE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE, 
              main = paste("Differentially expressed genes between \n all pairwise comparisons \n n = ", text1 , sep =))


plot_q

# Tip: If euclidean doesn't look good, try "correlation" other clustering distances
# Tip: If complete doesn't look good, try other methods
# Save the clustered heatmap plot as a .pdf:

#Plot Heatmap #1:
pdf("../03_output/HM1_clustered_genes.pdf", width = 6, height = 8)

plot_q

dev.off()
dev.off()

############### SPLITCONDITIONS: Split out heatmaps ################


# These are all the differentially expressed genes in all pair-wise comparisons:
summary(results_in_leaves_2)
summary(results_in_petioles_2)
summary(results_EV_2)
summary(results_HopBA1_2)
#Subset each results table for just the differentially expressed genes: Higher stringency
sign_leaf <- subset(results_in_leaves_2, padj < 0.001)
sign_petiole <- subset(results_in_petioles_2, padj < 0.001)
sign_EV <- subset(results_EV_2, padj < 0.001)
sign_Hop <- subset(results_HopBA1_2, padj < 0.001)

dim(sign_leaf)
dim(sign_petiole)
dim(sign_EV)
dim(sign_Hop)

#Determine how many genes were differentially expressed between HopBA1 v. EV in either leaf or petiole contexts, and merge them:
changing_genes <- rbind(sign_leaf, sign_petiole)

dim(changing_genes)
length(unique(rownames(changing_genes)))
# Get all r-stabilized log values of all the data:
rdl_all <- assay(rlog(CCdds, blind=FALSE))

#Subset just the changing genes (same as before just to get a smaller matrix):
changing_lrt_rdl <- subset(rdl_all, rownames(rdl_all) %in% rownames(changing_genes))
dim(changing_lrt_rdl)
head(changing_lrt_rdl)

## Heat Map #2

#Subset just the leaf data - genes changing in leaves - leaf samples only: 
changing_lrt_rdl_LEAF <- changing_lrt_rdl[,c(1,2,3,7,8,9)]

# Make sure it is in matrix form:
class(changing_lrt_rdl_LEAF)

#Subset the LEAF differential genes
changing_lrt_rdl_LEAF <- subset(changing_lrt_rdl_LEAF, rownames(changing_lrt_rdl_LEAF) %in% rownames(sign_leaf))
text2 = dim(changing_lrt_rdl_LEAF)[1]

# Draw Heatmap #2
# Scale by row, listed below as scale = "row". This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:
plot_r <- pheatmap(changing_lrt_rdl_LEAF, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=TRUE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE, 
              main = paste("Differentially expressed genes between HopBA1 and EV \n within leaves\n n = ", text2, sep = ""))

plot_r

# Save Heatmap #2 as a .pdf

pdf("../03_output/HM2_leaf_clustered_genes.pdf", width = 6, height = 8)

plot_r

dev.off()
dev.off()

## PETIOLE ###

#Subset just the petiole data:
changing_lrt_rdl_PETIOLE <- changing_lrt_rdl[,c(4,5,6,10,11,12)]
dim(changing_lrt_rdl_PETIOLE)
head(changing_lrt_rdl_PETIOLE)


# Make sure it is in matrix form:
class(changing_lrt_rdl_PETIOLE)

#Subset the LEAF data
changing_lrt_rdl_PETIOLE <- subset(changing_lrt_rdl_PETIOLE, rownames(changing_lrt_rdl_PETIOLE) %in% rownames(sign_petiole))
text3 = dim(changing_lrt_rdl_PETIOLE)[1]
text3
# Draw a heat map
# Scale by row, listed below as scale = "row". This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:
plot_s <- pheatmap(changing_lrt_rdl_PETIOLE, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=TRUE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE,
              main = paste("Differentially expressed genes between HopBA1 and EV \n within petioles\n n = ", text3, sep = ""))

plot_s


# Tip: If euclidean doesn't look good, try "manhattan" other clustering distances
# Tip: If complete doesn't look good, try other methods
# Save the clustered heatmap plot as a .pdf:

pdf("../03_output/HM3_petiole_clustered_genes.pdf", width = 6, height = 8)

plot_s

dev.off()
dev.off()


############## UPSET PLOT #################

# I want to find the intersection between different differential expression categories
dim(sign_leaf)
dim(sign_petiole)
head(sign_leaf)
head(sign_petiole)



# Get versions
sessionInfo()

# Please cite R in your publications when using it
citation()

