tnbcHRDscar R package manual
========================

-   [Introduction](#introduction)
    -   [Contact](#contact)
-   [Before to use it](#before-to-use-it)
    -   [Requirements](#requirements)
    -   [Installation](#installation)
    -   [Citation](#citation)
-   [Hands-on](#hands-on)
    -   [Input file example](#input-file-example)
    -   [Usage example](#usage-example)
-   [References](#references)

Introduction
============

R package to quantify specific allelic imbalances associated with homologous recombination deficiency in triple negative breast cancer (tnbc). Use as input a list of allele specific copy number segments per samples as those generated by the software ASCAT.

The specific allelic imbalances quantified by this package are:

1. LOH events longer than 10Mb but smaller than 30Mb
2. Large scale transitions (LSTs): Two consecutive allelic imbalances of at least 5Mb, with a separation between then no more than 2Mb
3. Telomeric allelic imbalances: Allelic imbalances longer than 1Mb that extend up to the telomere but do not cross the centromer.

The total number of the previous allelic imbalances are the resultant tnbcHRDscar value.

Contact
------------
- Fernando Perez-Villatoro (fernando.perez@helsinki.fi)


Before to use it
============

Requirements
-----------------
- Installation of R in any Operating system: Linux, OS X or Windows
- R recomended version 4.3.0
- R package `devtools`:
``` r
install.packages("devtools")
```

Installation
------------

The  `tnbcHRDscar` can be installed via R::`devtools` from github:

``` r
library(devtools)
install_github('farkkilab/tnbcHRDscar')
```

Citation
--------

Please include this citation to your references if you use this tool:

*Perez-Villatoro, F., Oikkonen, J., Casado, J., Chernenko, A., Gulhan, D. C., Tumiati, M., ... & Färkkilä, A. (2022). Optimized detection of homologous recombination deficiency improves the prediction of clinical outcomes in cancer. NPJ precision oncology, 6(1), 96.*
DOI: https://doi.org/10.1038/s41698-022-00339-8

Hands-on
=================

Input file example
-------------------
The input is a list of allele specific copy number (ASCN) segments generated by tools like [ASCAT](https://www.crick.ac.uk/research/labs/peter-van-loo/software). Example of input file:



``` r
a<-read.table("/examples/segments.txt", header=T)
head(a)
```

    ##         SampleID Chromosome Start_position End_position total_cn A_cn B_cn
    ## 1 SamplePatient1       chr1          14574       952448        5    0    5
    ## 2 SamplePatient1       chr1         953394      1259701        3    0    3
    ## 3 SamplePatient1       chr1        1278085      4551743        2    0    2
    ## 4 SamplePatient1       chr1        4551885     14124232        2    0    2
    ## 5 SamplePatient1       chr1       14161231     31062374        3    1    2
    ## 6 SamplePatient1       chr1       31074785     47428120        4    2    2

The input list should contain columns with the previous order. Each row represent a ASCN segment. It is not necessary to keep a proper row order.  

Column  description:

- `A_cn`: A allele copy number value

- `B_cn`: B allele  copy number value

- `total_cn`: Sum of A and B values

**Note: For each segment, the package will re-order the A and B copy number values, considering as A the one with highest ASCN.*

Usage example
-------------
Get the number of allelic imbalances associated with HRD in triple negative breast cancer with the function `get.tnbcHRDscars`:
``` r
library("tnbcHRDscar")
a <-read.table("/examples/segments.txt", header=T)
get.tnbcHRDscars(a, chrominfo ="grch38")
```

After running it will produce the next output:

    ##           nLOH LSTs nTAIs tnbcHRDscar
    ## SamplePatient1  25  35    33        93


Ruing parameters:
- `a` -- Input dataframe with segments  
- `chrominfo` -- the reference genome used, `grch38` or `grch37` (default: `grch38`)  

References
==========

*Some of the functions used are adaptations of the package scarHRD (https://github.com/sztup/scarHRD).* 

