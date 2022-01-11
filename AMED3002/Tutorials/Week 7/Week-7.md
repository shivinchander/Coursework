Week 7
================
Shivin Chander
1/11/2022

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

## Install / load required R packages

``` r
library(tidyverse)
library(ggfortify)
library(RColorBrewer)

# BiocManager needed for the next packages
# if (!requireNamespace('BiocManager', quietly =
# TRUE)) install.packages('BiocManager') /// BiocManager::install('PACKAGENAME')

library(GEOquery)
library(org.Hs.eg.db)
library(limma)
library(edgeR)
library(Glimma)
library(Homo.sapiens)
```

## Q1: Read in the data.

``` r
sfiles <- getGEOSuppFiles("GSE126379")
fnames <- rownames(sfiles)

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

## Q2: Tidy up the dataset and set it up for downstream analysis

``` r
ENTREZID <- mapIds(Homo.sapiens, keys = as.character(Data$Gene.Name), column = c("ENTREZID"), 
    keytype = "SYMBOL", multiVals = "first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
genes <- data.frame(ENTREZID, SYMBOL = Data$Gene.Name)

chooseGenes <- which(!duplicated(genes$SYMBOL) & !duplicated(genes$ENTREZID))

genes <- genes[chooseGenes, ]
Data <- Data[chooseGenes, ]
```

``` r
# MANIPULATE GENE COUNTS

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

# LCL are labelled 1, B cell labelled 2. renaming the columns to reflect this
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
# make the Gene IDs the name of the rows
rownames(counts) <- genes$SYMBOL
```

``` r
# MANIPULATE SAMPLE ANNOTATIONS

str_split_n <- function(string, pattern, n) {
    out <- str_split(string, pattern, simplify = TRUE)
    apply(out, 1, `[`, i = n)
}

# using sample names to make factors for group and individual
group <- str_split_n(colnames(counts), "_", 2) %>% factor()
individual <- str_split_n(colnames(counts), "_", 1) %>% factor()

# put information into a data frame
samples <- data.frame(group, individual)
samples
```

    ##    group individual
    ## 1    LCL          A
    ## 2  CD19B          A
    ## 3    LCL          B
    ## 4  CD19B          B
    ## 5    LCL          C
    ## 6  CD19B          C
    ## 7    LCL          D
    ## 8  CD19B          D
    ## 9    LCL          E
    ## 10 CD19B          E

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
# # calculate normalization factors
# BLCL <- calcNormFactors(BLCL, method = "TMM")
# # make a copy of BLCLDGEList to use in extension tab.
# BLCL2 <- BLCL
# 
# # show the nomalisation factors. For this dataset the effect of TMM-normalisation
# # is mild, as evident in the magnitude of the scaling factors, which are all
# # relatively close to 1.
# 
# BLCL$samples$norm.factors
# 
# # keep genes with large enough counts to be considered useful.
# keep <- filterByExpr(BLCL)
# BLCL <- BLCL[keep, ]
# 
# # calculate the counts per million reads for each gene for all samples
# cpm <- cpm(BLCL)
# # calculate the log of these CPM values
# lcpm <- cpm(BLCL, log = TRUE)
```

## Q4: Explore data with PCA. Can you see any structure?

``` r
# # Perform PCA
# pca <- prcomp(t(lcpm), scale = TRUE)
# # Plot results
# autoplot(pca, data = samples, colour = "group", shape = "individual")
# 
# # dimension reduction called MDS.
# plotMDS(lcpm, labels = group)
# plotMDS(lcpm, labels = individual)
# 
# 
# # using Glimma
# glMDSPlot(lcpm, labels = paste(group, individual, sep = "_"), groups = BLCL$samples[, 
#     c("group", "individual")], launch = TRUE)
```
