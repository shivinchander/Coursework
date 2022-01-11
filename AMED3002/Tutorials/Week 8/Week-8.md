Week 8
================
Shivin Chander
4/28/2021

# 2.1 B cells and EBV infected B cells (LCLs)

The following data comes from a manuscript “Evidence from genome wide
association studies implicates reduced control of Epstein-Barr virus
infection in multiple sclerosis susceptibility”
(<https://www.ncbi.nlm.nih.gov/pubmed/31039804>). There were multiple
goals for the RNAseq component of project:

  - Identify the gene expression changes that occur when B cells are
    infected with EBV and transformed into LCLs.

  - Determine which pathways are over-represented in the differentially
    expressed genes.

  - Assess if there are more MS risk genes then expected by chance in
    the differentially expressed genes.

Use this dataset as an opportunity to explore the behaviour of a small
RNA-Seq experiment. In this lab we will start by just processing and
exploring the data.

IMPORTANT: The code from week 8 continues on from the code generated in
week 7. The new steps introduced this week start from step 5.

``` r
# if (!requireNamespace('BiocManager', quietly = TRUE)) 
#     install.packages('BiocManager')
```

``` r
library(tidyverse)
library(ggfortify)
library(RColorBrewer)


# If any packages of these required are not installed, they can be installed by
# first installing bioconductor: if (!requireNamespace('BiocManager', quietly =
# TRUE)) install.packages('BiocManager') then using
# BiocManager::install('PACKAGENAME')

# load packages to be used in the analysis
library(GEOquery)
library(org.Hs.eg.db)
library(limma)
library(edgeR)
library(Glimma)
library(Homo.sapiens)
library(gplots)
library(dplyr)
```

``` r
# The count data for this GEO submission is available as a supplementary file on
# GEO. You only want to run this once...

sfiles <- getGEOSuppFiles("GSE126379")
fnames <- rownames(sfiles)

# There is only one supplemental file
Data = read.delim(fnames[1], header = TRUE)
head(Data)
```

    ##   Chromosome Gene.Name Gene.ID Sum.of.Exon.Length.of.Gene
    ## 1      chr10      A1CF    A1CF                       9529
    ## 2      chr10     ABCC2   ABCC2                       5051
    ## 3      chr10      ABI1    ABI1                       3776
    ## 4      chr10    ABLIM1  ABLIM1                       8303
    ## 5      chr10    ACADSB  ACADSB                       5941
    ## 6      chr10     ACBD5   ACBD5                       4066
    ##   Reads.Counts..BLCL_A1. RPKM..BLCL_A1. Reads.Counts..BLCL_A2. RPKM..BLCL_A2.
    ## 1                      2      0.0113048                      1    0.009967274
    ## 2                     21      0.2239352                     42    0.789760929
    ## 3                   1975     28.1718559                   1436   36.119869850
    ## 4                    786      5.0987985                   2619   29.958784090
    ## 5                    756      6.8539761                    543    8.680884936
    ## 6                    535      7.0870705                    421    9.834186588
    ##   Reads.Counts..BLCL_B1. RPKM..BLCL_B1. Reads.Counts..BLCL_B2. RPKM..BLCL_B2.
    ## 1                      1     0.01277064                      4     0.02185651
    ## 2                     10     0.24092549                     89     0.91744656
    ## 3                    563    18.14414625                   2307    31.81147083
    ## 4                    220     3.22439153                   4464    27.99350563
    ## 5                    210     4.30149944                    877     7.68613807
    ## 6                    200     5.98580756                    674     8.63099098
    ##   Reads.Counts..BLCL_C1. RPKM..BLCL_C1. Reads.Counts..BLCL_C2. RPKM..BLCL_C2.
    ## 1                      0      0.0000000                     11       0.122567
    ## 2                     13      0.2041644                     63       1.324315
    ## 3                    695     14.6004734                   1235      34.726663
    ## 4                    400      3.8215458                   2600      33.248066
    ## 5                    264      3.5249949                    568      10.151188
    ## 6                    304      5.9308963                    492      12.847712
    ##   Reads.Counts..BLCL_D1. RPKM..BLCL_D1. Reads.Counts..BLCL_D2. RPKM..BLCL_D2.
    ## 1                      2     0.02147553                     10     0.04834247
    ## 2                     25     0.50643522                    163     1.48657342
    ## 3                    847    22.95158503                   2806    34.23196323
    ## 4                    626     7.71437158                   5804    32.20093570
    ## 5                    272     4.68457946                   1065     8.25783462
    ## 6                    250     6.29120581                   1139    12.90424131
    ##   Reads.Counts..BLCL_E1. RPKM..BLCL_E1. Reads.Counts..BLCL_E2. RPKM..BLCL_E2.
    ## 1                      1     0.01387007                      7     0.04338079
    ## 2                     54     1.41300065                    130     1.51989185
    ## 3                    986    34.51206142                   2471    38.64447949
    ## 4                    177     2.81750179                   4484    31.89167049
    ## 5                    302     6.71851610                   1039    10.32767820
    ## 6                    313    10.17426241                    932    13.53615256

``` r
# Next time you can just run... fnames <-
# 'GSE126379/GSE126379_BLCL_processed.txt.gz' Data =
# read.delim(fnames[1],header=TRUE) head(Data)
```

## Q2: Tidy up the dataset and set it up for downstream analysis

``` r
# MANIPULATE GENE ANNOTATIONS

# Genes have multiple 'identifiers' in different databases. Later on we will need
# identifiers from the ENTREZ database. To save us some pain later, lets extract
# these now using the select() function from the AnnotationDbi package. The
# ENTREZ identifiers can be found in the Homo.sapiens package.

ENTREZID <- mapIds(Homo.sapiens, keys = as.character(Data$Gene.Name), column = c("ENTREZID"), 
    keytype = "SYMBOL", multiVals = "first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes <- data.frame(ENTREZID, SYMBOL = Data$Gene.Name)

# In this file of counts, there are multiple instances where the same gene symbol
# or entrez id appears in multiple rows. This may be due to how it was initially
# generated but as we are not able to distinguish between them (and in most cases
# the count information is very similar for the duplicate rows, we will just ask
# for just one row to be kept for any given gene name)

chooseGenes <- which(!duplicated(genes$SYMBOL) & !duplicated(genes$ENTREZID))

genes <- genes[chooseGenes, ]
Data <- Data[chooseGenes, ]


# MANIPULATE GENE COUNTS

# We do not need the RPKMs just the read counts, so just pull out these
# columns
colnames(Data)
```

    ##  [1] "Chromosome"                 "Gene.Name"                 
    ##  [3] "Gene.ID"                    "Sum.of.Exon.Length.of.Gene"
    ##  [5] "Reads.Counts..BLCL_A1."     "RPKM..BLCL_A1."            
    ##  [7] "Reads.Counts..BLCL_A2."     "RPKM..BLCL_A2."            
    ##  [9] "Reads.Counts..BLCL_B1."     "RPKM..BLCL_B1."            
    ## [11] "Reads.Counts..BLCL_B2."     "RPKM..BLCL_B2."            
    ## [13] "Reads.Counts..BLCL_C1."     "RPKM..BLCL_C1."            
    ## [15] "Reads.Counts..BLCL_C2."     "RPKM..BLCL_C2."            
    ## [17] "Reads.Counts..BLCL_D1."     "RPKM..BLCL_D1."            
    ## [19] "Reads.Counts..BLCL_D2."     "RPKM..BLCL_D2."            
    ## [21] "Reads.Counts..BLCL_E1."     "RPKM..BLCL_E1."            
    ## [23] "Reads.Counts..BLCL_E2."     "RPKM..BLCL_E2."

``` r
counts <- Data[, grep("Counts", colnames(Data))]

# LCL are labelled 1, B cell labelled 2. rename the columns to reflect this
colnames(counts) <- sub("Reads.Counts..BLCL_", "", colnames(counts))
colnames(counts) <- sub("1", "LCL", colnames(counts))
colnames(counts) <- sub("2", "CD19B", colnames(counts))
colnames(counts) <- sub("CD19B.", "_CD19B", colnames(counts))
colnames(counts) <- sub("LCL.", "_LCL", colnames(counts))
colnames(counts)
```

    ##  [1] "A_LCL"   "A_CD19B" "B_LCL"   "B_CD19B" "C_LCL"   "C_CD19B" "D_LCL"  
    ##  [8] "D_CD19B" "E_LCL"   "E_CD19B"

``` r
# MAKE DGEList

# make DGEList object to connect gene counts with sample information
# BLCL <- DGEList(counts = counts, samples = samples, genes = genes)
# BLCL
```

## Q3: Further processing and normalisation of the data

``` r
# # COUNTS NORMALIZATION
# 
# # First calculate normalization factors
# BLCL <- calcNormFactors(BLCL, method = "TMM")
# 
# # make a copy of BLCLDGEList to use in extension tab.
# BLCL2 <- BLCL
# 
# # show the nomalisation factors. For this dataset the effect of TMM-normalisation
# # is mild, as evident in the magnitude of the scaling factors, which are all
# # relatively close to 1.
# 
# BLCL$samples$norm.factors
```

## Q4: Explore data with PCA. Can you see any structure?

``` r
# # Perform PCA
# pca <- prcomp(t(lcpm), scale = TRUE)
# # Plot results
# autoplot(pca, data = samples, colour = "group", shape = "individual")
```

``` r
# # Use another method for dimension reduction called MDS.
# plotMDS(lcpm, labels = group)
```

``` r
# plotMDS(lcpm, labels = individual)
```

``` r
#interactive plot of MDS using Glimma

# glMDSPlot(lcpm, labels = paste(group, individual, sep = "_"), groups = BLCL$samples[, 
#     c("group", "individual")], launch = TRUE)
```

## Q5: Construct a design matrix and perform variance stabilisation.

limma is a versatile package for DE analysis of microarray data.
However, the count based data from RNAseq has certain stastical
properties that are different to microarray based data, namely that in
RNAseq count data, the variance is not independent of the mean. The
authors of limma delevoped the voom function which accounts for this and
modifies the RNAseq count data for use with the limma package.

``` r
# We first need to set up our design matrix, where we tell limma that we are
# interested in fitting the gene expression data against group, but also to
# include the 'individual' term, as this will account for the samples being
# paired.

# design <- model.matrix(~individual + group)
# design
```

``` r
# v <- voom(BLCL, design, plot = TRUE)
# v
```

## Q6: Fit linear model

``` r
# Linear modelling in limma is carried out using the 'lmFit' and 'contrasts.fit'
# functions.  These functions fit a separate model to the expression values for
# each gene. We now tell limma to fit the model specified in our design matrix
# for the voom transformed data


# fit <- lmFit(v, design)
```

## Q7: Generate moderated test statistics

``` r
# efit <- eBayes(fit)
```

## Q8: Explore the results

``` r
# dt_fdr <- decideTests(efit, adjust.method = "fdr", p.value = 0.05)
# summary(dt_fdr)
```

``` r
# dt_bonf <- decideTests(efit, adjust.method = "bonferroni", p.value = 0.05)
# summary(dt_bonf)
```

``` r
# hist(efit$p.value[, "groupLCL"])
```

``` r
# volcanoplot(efit, coef = "groupLCL", highlight = 10, names = efit$genes$SYMBOL)
```

``` r
# DELCLvB <- topTable(efit, coef = "groupLCL", adjust.method = "fdr", p.value = 0.05, 
#     number = Inf)
```

``` r
# write.csv(DELCLvB, file = "DELCLvB.csv")
```

``` r
# useCoef <- which(colnames(efit) == "groupLCL")
# plotMD(efit, column = useCoef, status = dt_fdr[, useCoef], main = colnames(efit)[useCoef], 
#     )
```

``` r
# glMDPlot(efit, coef = useCoef, status = dt_fdr[, useCoef], main = colnames(efit)[useCoef], 
#     side.main = "ENTREZID", counts = lcpm, groups = group, launch = TRUE)
```

## Q9: Create a gene expression heatmap

``` r
# LCLvB.topgenes <- DELCLvB$SYMBOL[1:40]
# i <- which(v$genes$SYMBOL %in% LCLvB.topgenes)
# coolmap(lcpm[i, ], margin = c(8, 6), cexRow = 0.7, cexCol = 0.7)
# 
# 
# pdf("heatmap.pdf", width = 12, height = 16)
# coolmap(lcpm[i, ], margin = c(8, 6))
# dev.off()
```

## Q10: Draw a plot of a gene of interest from the data

``` r
# CD40 <- lcpm["CD40", ]
# df <- data.frame(CD40, group, individual)
# ggplot(df, aes(x = group, y = CD40, colour = individual)) + geom_point() + geom_line(aes(group = individual)) + 
#     theme_classic()
```
