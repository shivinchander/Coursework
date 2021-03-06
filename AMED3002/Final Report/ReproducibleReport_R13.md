Group Project 2
================
Remote Group 13
11/06/2021

``` r
library(ggplot2)
library(limma)
library(ggfortify)
library(tidyverse)
library(GEOquery)
library(MASS)
library(naniar)
library(dendextend)
library(org.Hs.eg.db)
library(edgeR)
library(Glimma)
library(Homo.sapiens)
library(gplots)
library(dplyr)
library(AnnotationDbi)
library(RColorBrewer)
library(e1071)
```

# Summary

## Aim

The aim of this study conducted by Wang and colleagues was to use both
transcriptomic databases and SVM-based pattern recognition to select
biomarkers that could be used in predicting and preventing gestational
diabetes mellitus (GDM) as well as determine links, if any, to other
correlating diseases such as Type 2 Diabetes and immune disorders.

## Methods

The study first screened 63 samples (32 GDM patient samples, 31 control
samples) from a dataset collected from GEO (GSE70493) for GDM-specific
biomarkers. The study then utilized the edgeR package to screen out the
differentially expressed genes between the two groups. The paper then
proceeded to use the Poisson regression model and Empirical Bayes method
to account for biological variability and improve the overall
reliability of future inferences.

Hierarchical clustering was then performed to create a heatmap
expressing the gene expression data analysis, and a volcano plot was
created in order to observe the significance of the microarrays gene
selection criteria. Through the use of databases such as DAVID (database
for annotation, visualization, and integration discovery) and KEGG
(Kyoto encyclopedia of genes and genomes), a gene regulatory network was
then constructed in order to select GDM biomarker candidates (utilizing
prior literature sourced from PubMed).

Once biomarker candidates were located, all were subjected to the
implementation of a SVM (support vector machine) to determine the
prediction capacity of the model. Finally, a 10-fold cross validation
was performed to determine the validity of the study’s findings and
observe the error rate.

## Results

A total of 6 genes were found to serve as good predictors of GDM and
thus have the potential to serve as biomarkers for GDM prediction and
diagnosis. The biomarkers found were all derived from the HLA
superfamily - HLA-DQA1, HLA-DQA2, HLA-DRA, HLA-DRB1, HLA-DRB4 and
HLA-DRB5.

# Experimental Design

The aforementioned methodologies of the paper can be split into 4
components:  
1\. Obtaining of gene expression profiles  
2\. Analysis of genes to narrow down potential biomarkers  
3\. Obtaining of potential GDM biomarkers  
4\. Creation of prediction models for validity

## Workflow

``` r
knitr::include_graphics("Capture.png")
```

<img src="Capture.png" width="1368" />

# Reproduced Report

## Load Data

``` r
gset <- getGEO("GSE70493",getGPL=TRUE, destdir = ".")[[1]]
```

## Rename Pheno

``` r
exprs <- exprs(gset)
geneAnno <- fData(gset)
pheno <- pData(gset)
view(pheno)
colnames(pheno) <- gsub(" ", "_", colnames(pheno))
colnames(pheno) <- sub(":ch1", "", colnames(pheno))

x <- log2(exprs)
```

## Remove Missing Data

``` r
Data <- na.omit(gset)
```

This was not explicitly stated that it was done in the current paper. We
decided to remove missing data anyway in order to get more accurate
results/findings to compare to the original papers own findings.

## Creation of Design Matrix and Consequent Fits

``` r
design <- model.matrix(~pheno$gestational_diabetes, pheno)

fit = lmFit(x, design)
efit = eBayes(fit, trend = TRUE)

transpose <- t(x)
dt_fdr <- decideTests(efit, method = "hierarchical", p.value = 0.05)
summary(dt_fdr)
```

    ##        (Intercept) pheno$gestational_diabetesYes
    ## Down             0                          1262
    ## NotSig           0                         68035
    ## Up           70523                          1226

This was an interesting finding within our reproduced report as within
the paper, they stated that they found 156 downregulated and 33
upregulated probes, making it a total of 189 differentially expressed
probes. In our report, we instead found 1,257 downregulated and 1,218
upregulated probes - obviously very different from the original.

Unfortunately, the methods section was not explicit enough to state how
they got these probes narrowed down and thus we could not replicate it
to the full extent. This, however, is an improvement on prior attempts
made by members of the group, as previously we had found that none of
the probes were significant, which would have been less helpful then our
current findings.

It is important to note that the paper also made use of both the Poisson
regression model and Empirical Bayes method. Whilst we did use Empirical
Bayes Method, as can be seen in the above coding, we decided not to use
Poisson regression model as there was no impact of applying it to the
data. This was confirmed by a demonstrator during a 1-on-1 session.

## Hierarchical Clustering and Heatmap

``` r
## Hierarchical Clustering

library(pheatmap)
sampleInfo <- pData(gset)
gene_matrix <- cor(exprs(gset),use="c")

pheatmap(gene_matrix,
     labels_row = geneAnno$probeset_id,
     scale="row")
```

![](ReproducibleReport_R13_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

The study created its heatmap from utilising a dendrogram created from
hierarchical clustering. This was done by isolating the differentially
expressed genes that had a p-value of less than 0.05 using edgeR.
Unfortunately, their methodologies on how this was performed (especially
with a dataset with 70,000 probes) was unavailable.

After multiple attempts of trying to minimise the list of genes utilised
and expressed on the heatmap, we were unable to reproduce this exactly,
finding that methodologies which could do this successfully online
utilised prior versions of rStudio which we did not have access to.
Furthermore, it was not explicitly stated which form of hierarchical
clustering was used within the data analysis.

Instead, we chose to create a heatmap, showing a different part of the
data available. The above heatmap shows a comparison of the probeset IDs
and samples, with the colour scheme representing the expression levels
of the differentially expressed genes, as found on the design matrix
created. Orange represents high expression levels and blue represents
low expression levels.

## Volcano Plot

``` r
-log10(0.05)
```

    ## [1] 1.30103

``` r
top100DE <- rownames(dt_fdr)[1:100]
top19DE <- rownames(dt_fdr)[1:19]
columnColour = pheno$gestational_diabetes
columnColour[pheno$gestational_diabetes == "Yes"] = "red"
columnColour[pheno$gestational_diabetes == "No"] = "white"
useColumns <- !is.na(columnColour)
```

## DAVID and KEGG

The manuscript, after these initial differential analyses, performed
both an enrichment analysis via the DAVID tool and constructed a KEGG
pathway network in order to create a gene regulatory network - in other
words, create gene symbols to replace the gene IDs provided in the
initial dataset.

For our reproduced report, we did not make use of the DAVID database as
upon inspection it required all 70,000 probes to be manually typed in to
observe the gene IDs. Furthermore, the KEGG pathway could also not be
conducted due to the paper not being explicit in its use and
contribution to the final results.

## SVM

``` r
tab <- topTable(efit, coef = 2, adjust = "fdr", n = 10)
labCol <- c("Yes", "No")
top10 <- topTable(efit)

X <- t(exprs[rownames(top10), ])  # removed colnames[design] // wasnt needed

# Select outcome variable.
Y <- as.factor(pheno[, "gestational_diabetes"]) # as above, removed colnames[design] from left of comma

# test <- factor(rep(NA, n), levels = c("No", "Yes")) # dont need :(

# create new dataset.
gset1 <- data.frame(Y, X)

fit <- svm(Y ~ ., gset1)

predSVM <- predict(fit, gset1)  # Here we have predicted the classes of our original data.

mean(predSVM != gset1$Y) # value of 0.1269841
```

    ## [1] 0.1746032

The current paper utilised a support vector machine (SVM) to predict the
potential application of the 6 found biomarkers in GDM clinical
diagnosis. Through the creation of a test set and training set
containing these particular data points, the study was able to obtain an
error rate of 29.03%. In our reproduction of this prediction model, we
got an error rate of 17.4%, which gives the impression that the dataset
is not as accurate as what would typically be ideal within the
scientific community. Despite this, it was still a better error rate
than the original report, which leads to questions regarding the
methodology within the paper and whether it misled our interpretation.

## randomForest

``` r
library(randomForest)
fit <- randomForest(Y ~ ., gset1)
predrandomForest <- predict(fit, gset1)  
predrandomForest <-factor(predrandomForest, levels = c("No", "Yes"))
predrandomForest
```

    ## GSM1784987 GSM1784988 GSM1784989 GSM1784990 GSM1784991 GSM1784993 GSM1784994 
    ##        Yes        Yes        Yes         No        Yes         No        Yes 
    ## GSM1784995 GSM1784996 GSM1784997 GSM1784998 GSM1784999 GSM1785000 GSM1785001 
    ##        Yes        Yes        Yes         No        Yes        Yes         No 
    ## GSM1785003 GSM1785004 GSM1785005 GSM1785006 GSM1785007 GSM1785008 GSM1785009 
    ##         No         No         No         No        Yes        Yes        Yes 
    ## GSM1785011 GSM1785012 GSM1785013 GSM1785014 GSM1785015 GSM1785016 GSM1785017 
    ##        Yes        Yes         No         No        Yes        Yes         No 
    ## GSM1785018 GSM1785020 GSM1785021 GSM1785022 GSM1785023 GSM1785024 GSM1785025 
    ##        Yes         No         No         No        Yes         No         No 
    ## GSM1785027 GSM1785028 GSM1785029 GSM1785030 GSM1785031 GSM1785032 GSM1785033 
    ##         No         No         No        Yes         No        Yes         No 
    ## GSM1785035 GSM1785036 GSM1785037 GSM1785038 GSM1785039 GSM1785040 GSM1785041 
    ##         No         No         No         No        Yes        Yes        Yes 
    ## GSM1785042 GSM1785044 GSM1785045 GSM1785046 GSM1785047 GSM1785048 GSM1785049 
    ##        Yes         No        Yes        Yes        Yes        Yes         No 
    ## GSM1785050 GSM1785052 GSM1785053 GSM1785054 GSM1785055 GSM1785056 GSM1785057 
    ##         No        Yes         No         No        Yes         No        Yes 
    ## Levels: No Yes

``` r
gset1$Y
```

    ##  [1] Yes Yes Yes No  Yes No  Yes Yes Yes Yes No  Yes Yes No  No  No  No  No  Yes
    ## [20] Yes Yes Yes Yes No  No  Yes Yes No  Yes No  No  No  Yes No  No  No  No  No 
    ## [39] Yes No  Yes No  No  No  No  No  Yes Yes Yes Yes No  Yes Yes Yes Yes No  No 
    ## [58] Yes No  No  Yes No  Yes
    ## Levels: No Yes

``` r
table(predrandomForest)
```

    ## predrandomForest
    ##  No Yes 
    ##  31  32

``` r
table(gset1$Y)
```

    ## 
    ##  No Yes 
    ##  31  32

``` r
mean(predrandomForest != gset1$Y)
```

    ## [1] 0

``` r
## Gave us zero so we do this...

set.seed(51773)
n <- nrow(gset1)
trainTest <- rep(c("train", "test"), 100)
trainTest <- trainTest[1:n]
trainTestSample <- sample(trainTest, n, replace = FALSE)
dataTrain <- gset1[trainTestSample == "train", ]
dataTest <- gset1[trainTestSample == "test", ]

fitforest <- randomForest(Y ~ ., dataTrain)
predrandomForest <- predict(fitforest, dataTest)
mean(predrandomForest != dataTest$Y)
```

    ## [1] 0.483871

To test the data provided using another prediction model, we decided to
use randomForest. RandomForest can process a mixture of numerical and
categorical features and does not require the data to be scaled, which
is to our advantage since we had over 70,000 probes to process. Our
findings indicated an error rate of 48.3%, which was greater than that
found using the SVM model. Due to the high rate of error, this shows
that the paper may have been misleading in its findings, this suggestion
supported by the issues we found throughout the methodology component.

# Appropriateness

## Methods

The study was able to successfully determine 6 biomarkers that could be
used clinically in the diagnosis of GDM, thus achieving their goal, thus
indicating that their methodology served its purpose and was correct in
order to give the provided outcome. The methodology of utilising the
DAVID database and KEGG pathway tools was found to be inappropriate as
these programs were found to be inaccessible and difficult to use,
especially when there was a requirement to manually insert the 70,000
probes - which in our time frame, was an unreachable goal.

Furthermore, there are doubts as to whether the methods have sufficient
explanation in order to be reproducible for peer-reviews. So despite
their benefit and success within this paper, it cannot be stated with
confidence that the methods were appropriate due to the lack of
reproducibility.

## Analysis

The analytical portion of the paper, particularly the SVM prediction
model and 10-fold cross validation were for the most-part appropriate as
they allowed for the determination of the validity of the model created
within the study. For the prediction model, SVM is not the most
appropriate choice for analysis due to the large scale of the data.
Prediction models such as randomForest are better equipped for such
datasets and thus would have potentially given a more accurate error
rate and overall perspective of the predictive model created for the six
biomarkers.

In regards to the 10-fold cross validation, we believe that this was an
appropriate analytical method, as it provides the study with greater
accuracy as to the fitting procedure as well as reduces both the
variance and bias. Therefore, deciding to use a 10-fold cross validation
over a 5-fold is ideal in this scenario.

# Successes and Difficulties

## Successes

Despite our difficulties in extrapolating and reproducing this paper, we
did have some successes. Our group was able to understand the
methodologies used, despite not being able to reproduce them, and follow
the line of reasoning as to why the paper made these decisions.
Furthermore, we were able to reproduce, though not accurately, the SVM
prediction model as well as implement our own randomForest classifier to
compare the error rates. Additionally, we were able to utilise the
dataset and formulate our own heatmaps and volcano plots. These may have
not been reproduced from the paper, but it was able to show that we
understood the data well enough to extract different variables from it
to analyse.

## Difficulties

As evident from the above reproduced report, we encountered many
difficulties with the paper. One of our greatest weaknesses in
reproducing the report was some of the databases and programs utilised
by the paper, such as DAVID and KEGGS, which were difficult and too
time-consuming to utilise as intended by the paper. In turn, this meant
that we were unable to create gene symbols within our reproduced report,
which led to many difficulties later on in the paper. Additionally, we
had difficulty in reproducing many components of the report due to the
ambiguous nature of the writing. There was a lack of sufficient detail
that would typically be required in order to reproduce someones findings
in the case of peer-reviews. This may have been due to the fact that
this dataset was not created by the authors of the paper we looked at,
rather they used another papers dataset and created their own
interpretation of the findings through analytical methods on rStudio.
This could have led to a mis-translation within the methods section,
thus making it difficult to reproduce.

overall, whilst a difficult paper to reproduce, we have learnt a lot
about the difficulties of data analysis and interpretation, especially
with such large datasets. It can be difficult to manage and organise, as
well as interpret. This project has given all team members a great
appreciation for the complexities of such work.

## Improvements to Workflow

Potential improvements to the workflow would include the use of other
programs within rStudio rather than DAVID and KEGG to avoid unneccessary
complexities within the methodology of extracting information from the
provided dataset. Whether this be through re-working the existing
dataset to create subsets or creating a clear pathway to extracting the
necessary gene symbols in order for the subsequent coding and findings
to be parallel to that of the paper.

# References

Wang, Y., Wang, Z., & Zhang, H. (2018). Identification of diagnostic
biomarker in patients with gestational diabetes mellitus based on
transcriptome‐wide gene expression and pattern recognition. Journal Of
Cellular Biochemistry, 120(2), 1503-1510.
<https://doi.org/10.1002/jcb.27279>

``` r
> sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] randomForest_4.6-14                     pheatmap_1.0.12                         e1071_1.7-4                             RColorBrewer_1.1-2                     
 [5] gplots_3.1.1                            Homo.sapiens_1.3.1                      TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 GO.db_3.12.1                           
 [9] OrganismDbi_1.32.0                      GenomicFeatures_1.42.3                  GenomicRanges_1.42.0                    GenomeInfoDb_1.26.7                    
[13] Glimma_2.0.0                            edgeR_3.32.1                            org.Hs.eg.db_3.12.0                     AnnotationDbi_1.52.0                   
[17] IRanges_2.24.1                          S4Vectors_0.28.1                        dendextend_1.14.0                       naniar_0.6.0                           
[21] MASS_7.3-53.1                           GEOquery_2.58.0                         Biobase_2.50.0                          BiocGenerics_0.36.0                    
[25] forcats_0.5.1                           stringr_1.4.0                           dplyr_1.0.5                             purrr_0.3.4                            
[29] readr_1.4.0                             tidyr_1.1.3                             tibble_3.1.0                            tidyverse_1.3.0                        
[33] ggfortify_0.4.11                        limma_3.46.0                            ggplot2_3.3.3                          

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                backports_1.2.1             BiocFileCache_1.14.0        plyr_1.8.6                  lazyeval_0.2.2             
  [6] splines_4.0.2               BiocParallel_1.24.1         digest_0.6.27               htmltools_0.5.1.1           viridis_0.5.1              
 [11] fansi_0.4.2                 magrittr_2.0.1              memoise_2.0.0               Biostrings_2.58.0           annotate_1.68.0            
 [16] modelr_0.1.8                matrixStats_0.58.0          askpass_1.1                 rmdformats_1.0.1            prettyunits_1.1.1          
 [21] colorspace_2.0-0            blob_1.2.1                  rvest_1.0.0                 rappdirs_0.3.3              haven_2.3.1                
 [26] xfun_0.22                   crayon_1.4.1                RCurl_1.98-1.3              jsonlite_1.7.2              graph_1.68.0               
 [31] genefilter_1.72.1           survival_3.2-7              glue_1.4.2                  gtable_0.3.0                zlibbioc_1.36.0            
 [36] XVector_0.30.0              DelayedArray_0.16.3         DEoptimR_1.0-8              scales_1.1.1                DBI_1.1.1                  
 [41] Rcpp_1.0.6                  viridisLite_0.3.0           xtable_1.8-4                progress_1.2.2              bit_4.0.4                  
 [46] DT_0.17                     htmlwidgets_1.5.3           httr_1.4.2                  ellipsis_0.3.1              pkgconfig_2.0.3            
 [51] XML_3.99-0.5                dbplyr_2.1.0                locfit_1.5-9.4              utf8_1.1.4                  tidyselect_1.1.0           
 [56] rlang_0.4.10                reshape2_1.4.4              munsell_0.5.0               cellranger_1.1.0            tools_4.0.2                
 [61] cachem_1.0.4                cli_2.3.1                   generics_0.1.0              RSQLite_2.2.3               broom_0.7.5                
 [66] evaluate_0.14               fastmap_1.1.0               cvTools_0.3.2               yaml_2.2.1                  knitr_1.31                 
 [71] bit64_4.0.5                 fs_1.5.0                    robustbase_0.93-7           caTools_1.18.1              visdat_0.5.3               
 [76] RBGL_1.66.0                 xml2_1.3.2                  biomaRt_2.46.3              compiler_4.0.2              rstudioapi_0.13            
 [81] plotly_4.9.3                curl_4.3                    reprex_1.0.0                geneplotter_1.68.0          stringi_1.5.3              
 [86] lattice_0.20-41             Matrix_1.2-18               vctrs_0.3.6                 pillar_1.5.1                lifecycle_1.0.0            
 [91] BiocManager_1.30.12         data.table_1.14.0           bitops_1.0-6                rtracklayer_1.49.5          R6_2.5.0                   
 [96] bookdown_0.21               KernSmooth_2.23-18          gridExtra_2.3               gtools_3.8.2                assertthat_0.2.1           
[101] SummarizedExperiment_1.20.0 openssl_1.4.3               DESeq2_1.30.1               withr_2.4.1                 GenomicAlignments_1.26.0   
[106] Rsamtools_2.6.0             GenomeInfoDbData_1.2.4      hms_1.0.0                   grid_4.0.2                  class_7.3-18               
[111] rmarkdown_2.7               MatrixGenerics_1.2.1        lubridate_1.7.10  
> 
```

# INDIVIDUAL - 10-Fold Cross Validation (Method Utilised in Paper) - Shivin Chander

The paper uses a 10-Fold Cross Validation. Using certain cross
validation techniques such as single hold out uses a certain percent of
the data as training and the rest as the test set. For example, using
90% for training and 10% for testing, the test set is small and can
therefore have large amounts of variation in the estimates for different
samples of the data due to only fitting the data one time. 10-fold
validation helps reduce this variance by training on a random sample of
90% of the data and testing on the remaining 10% 10 times and averaging
out the results. Additionally, by using k-fold cross validation,
different partitioning of the dataset can be taken to form k sub-sets
and then averaged. This method can further increase the robustness of
the model. It should be noted that all steps of the model fitting must
be performed independently for each fold to reduce bias.

``` r
set.seed(666)
n <- nrow(gset1)  # The number of observations
nfolds <- 10

# Create a vector with the same amount of fold labels
fold <- rep(1:nfolds, 100)
fold <- fold[1:n]

# Reorder these to avoid systematic ordering bias
foldSample <- sample(fold, n, replace = FALSE)
```

``` r
# If response isn't a factor rep(NA,n)
predSVM <- factor(rep(NA, n), levels = c("Yes", "No"))
```

``` r
for (i in 1:nfolds) {
    
    dataTest <- gset1[foldSample == i, ]
    dataTrain <- gset1[foldSample != i, ]
    
    fit <- svm(Y ~ ., dataTrain)
    predSVM[foldSample == i] <- predict(fit, dataTest)
    
}
```

``` r
mean(predSVM != gset1$Y)
```

    ## [1] 0.2380952
