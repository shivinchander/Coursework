Week 3
================
Shivin Chander
3/17/2021

# 3.1 Respiratory illness

This data contains the respiratory status of patients recruited for a
randomised clinical multicenter trial. In each of two centres, eligible
patients were randomly assigned to active treatment or placebo. During
the treatment, the respiratory status (categorised poor or good) was
determined at each of four, monthly visits. The trial recruited 111
participants (54 in the active group, 57 in the placebo group) and there
were no missing data for either the responses or the covariates. The
question of interest is to assess whether the treatment is effective and
to estimate its effect.

### Q1: Read in the data and summarise it into a contingency table.

``` r
respiratory <- read.delim("https://wimr-genomics.vip.sydney.edu.au/AMED3002/data/respiratory.txt", 
    sep = "\t")

tab <- table(respiratory$treatment, respiratory$status)
tab
```

    ##            
    ##             good poor
    ##   placebo    127  158
    ##   treatment  172   98

### Q2: We would like to test if there is any evidence that receiving the treatment altered your chance of having a good respiratory status? Rephrase this question into a null and alternate hypothesis that is consistent with a chi-square test.

N0: There is no difference between the people treated and those that
were not.

Na: There is a difference between the people treated and those that were
not.

### Q3: Perform a chi-square test

``` r
chisq.test(tab)
```

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  tab
    ## X-squared = 19.682, df = 1, p-value = 9.148e-06

### Q4: Check the assumptions for a chi-square test

``` r
test = chisq.test(tab)
test$expected >= 5
```

    ##            
    ##             good poor
    ##   placebo   TRUE TRUE
    ##   treatment TRUE TRUE

``` r
str(test)
```

    ## List of 9
    ##  $ statistic: Named num 19.7
    ##   ..- attr(*, "names")= chr "X-squared"
    ##  $ parameter: Named int 1
    ##   ..- attr(*, "names")= chr "df"
    ##  $ p.value  : num 9.15e-06
    ##  $ method   : chr "Pearson's Chi-squared test with Yates' continuity correction"
    ##  $ data.name: chr "tab"
    ##  $ observed : 'table' int [1:2, 1:2] 127 172 158 98
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:2] "placebo" "treatment"
    ##   .. ..$ : chr [1:2] "good" "poor"
    ##  $ expected : num [1:2, 1:2] 154 145 131 125
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:2] "placebo" "treatment"
    ##   .. ..$ : chr [1:2] "good" "poor"
    ##  $ residuals: 'table' num [1:2, 1:2] -2.14 2.2 2.31 -2.38
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:2] "placebo" "treatment"
    ##   .. ..$ : chr [1:2] "good" "poor"
    ##  $ stdres   : 'table' num [1:2, 1:2] -4.52 4.52 4.52 -4.52
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:2] "placebo" "treatment"
    ##   .. ..$ : chr [1:2] "good" "poor"
    ##  - attr(*, "class")= chr "htest"

### Q5: What is your conclusion for this test?

The p-value calculated was 9.148e-06, therefore, the treatment did have
a positive association with producing a good result.

### Q6: Interpret the relationship further by calculating a relative risk or odds ratio. Are both appropriate in this case?

``` r
OR <- (tab[1, 1] * tab[2, 2])/(tab[2, 1] * tab[1, 2])
OR
```

    ## [1] 0.4579776

# 3.2 Random data (Extension)

Lets perform a test of independence with randomly generated data

### E1: Create two vectors of random data.

``` r
## Set seed for reproducible results
set.seed(51773)

## Generate two random vectors of size n. Can be done with the 'sample'
## function.
n = 50
AB = sample(c("A", "B"), n, replace = TRUE)
CD = sample(c("C", "D"), n, replace = TRUE)

head(AB)
```

    ## [1] "B" "B" "B" "B" "A" "A"

``` r
head(CD)
```

    ## [1] "C" "D" "D" "C" "C" "C"

### E2: Create a contingency table to view the relationship between AB and CD

``` r
tabABCD = table(AB, CD)
tabABCD
```

    ##    CD
    ## AB   C  D
    ##   A 18 10
    ##   B 11 11

### E3: Test for independence between AB and CD

``` r
chisq.test(tabABCD)
```

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  tabABCD
    ## X-squared = 0.529, df = 1, p-value = 0.467

### E4: Check assumptions of test

``` r
chisq.test(tabABCD)$expected
```

    ##    CD
    ## AB      C     D
    ##   A 16.24 11.76
    ##   B 12.76  9.24

### E5: Conclusion.

The p-value is greater than 0.05 therefore we fail to reject the null
hypothesis and we conclude that the result is statistically
nonsignificant. There is no evidence that AB and CD are not independent.
