# CORE

This is the code associated with the paper : 

Perturbations of dietary, gut mucosal, and microbial-metabolic signatures persist during remission in Crohn’s disease despite effective immune suppression


Background and Aims: We aim to identify dietary and gut signals during Crohn’s disease (CD) remission that differ from the healthy state to facilitate improved disease control and clearance.
Methods: We analyze diet, ileal transcriptomics, microbiomics, and metabolomics to compare patients with CD in remission to those with active CD, using non-IBD controls as the reference for healthy signals.
Results: 191 subjects were included; 77 CD patients in remission, 37 with active CD, and 77 non-IBD controls. Ileal transcriptomics revealed a significant decrease in genes and pathways associated with adaptive T-cells and innate granulocytes during remission, even lower than observed in non-IBD controls.  Despite the observed immune suppression, patients in remission showed an increase in the expression of epithelial antimicrobial pathways and related genes, including DUOX2, along with a rise in genes associated with goblet cells and mucin glycosylation. These signals during remission coincided with a persistent pathogenic gut microbial composition, metabolic alterations, and less healthy dietary habits, which were characterized by a higher intake of ultra-processed foods and lower consumption of fiber, folate, vitamin C, and vegetables. Greater exposure to ultra-processed foods was significantly associated with more dysbiotic gut signals and negatively correlated with genes enriched for mucin glycosylation, which is essential for maintaining gut barrier homeostasis.
Conclusion: Perturbations in dietary, ileal, microbial, and metabolic signatures persist during CD remission despite effective immune treatments. These findings likely underlie the remitting-relapsing nature of CD and suggest that interventions targeting diet, epithelial health, microbial, and metabolic functions may promote deeper, longer-lasting remission states.
 


The functions used through the code should all be available under the funcs/ directory (with some additional functions that were not eventually used here).  
Example data files are under the data/ subdirectory within the appropriate directory.


### Raw data is publicly available at the following:

The 16S amplicon sequencing dataset generated in this study has been deposited in BioProject under accession code: [PRJNA1288731](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1288731).


## R session info:
```
setting  value
 version  R version 4.1.2 (2021-11-01)
 os       Ubuntu 20.04.3 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Etc/UTC
 date     2025-07-15
 rstudio  2025.05.1+513 Mariposa Orchid (server)
 pandoc   3.4 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
 quarto   1.6.42 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto

─ Packages ────────────────────────────────────────────────────────────────────────────
 package              * version   date (UTC) lib source
 annotate               1.72.0    2021-10-26 [1] Bioconductor
 AnnotationDbi          1.56.2    2021-11-09 [1] Bioconductor
 backports              1.5.0     2024-05-23 [1] RSPM (R 4.1.0)
 base64enc              0.1-3     2015-07-28 [1] RSPM (R 4.1.0)
 biglm                  0.9-3     2024-06-12 [1] RSPM (R 4.1.0)
 Biobase              * 2.54.0    2021-10-26 [1] Bioconductor
 BiocGenerics         * 0.40.0    2021-10-26 [1] Bioconductor
 BiocParallel           1.28.3    2021-12-09 [1] Bioconductor
 Biostrings             2.62.0    2021-10-26 [1] Bioconductor
 bit                    4.6.0     2025-03-06 [1] RSPM (R 4.1.0)
 bit64                  4.6.0-1   2025-01-16 [1] RSPM (R 4.1.0)
 bitops                 1.0-9     2024-10-03 [1] RSPM (R 4.1.0)
 blob                   1.2.4     2023-03-17 [1] RSPM (R 4.1.0)
 cachem                 1.1.0     2024-05-16 [1] RSPM (R 4.1.0)
 caTools                1.18.3    2024-09-04 [1] RSPM (R 4.1.0)
 checkmate              2.3.2     2024-07-29 [1] RSPM (R 4.1.0)
 cli                    3.6.4     2025-02-13 [1] RSPM (R 4.1.0)
 cluster                2.1.2     2021-04-17 [2] CRAN (R 4.1.2)
 codetools              0.2-18    2020-11-04 [2] CRAN (R 4.1.2)
 colorspace             2.1-1     2024-07-26 [1] RSPM (R 4.1.0)
 crayon                 1.5.3     2024-06-20 [1] RSPM (R 4.1.0)
 curl                   6.2.1     2025-02-19 [1] RSPM (R 4.1.0)
 data.table             1.17.0    2025-02-22 [1] RSPM (R 4.1.0)
 DBI                    1.2.3     2024-06-02 [1] RSPM (R 4.1.0)
 DelayedArray           0.20.0    2021-10-26 [1] Bioconductor
 DEoptimR               1.1-3-1   2024-11-23 [1] RSPM (R 4.1.0)
 DESeq2               * 1.34.0    2021-10-26 [1] Bioconductor
 devtools             * 2.4.5     2022-10-11 [1] RSPM (R 4.1.0)
 digest                 0.6.37    2024-08-19 [1] RSPM (R 4.1.0)
 doParallel             1.0.17    2022-02-07 [1] RSPM (R 4.1.0)
 dplyr                  1.1.4     2023-11-17 [1] RSPM (R 4.1.0)
 dynamicTreeCut       * 1.63-1    2016-03-11 [1] RSPM (R 4.1.0)
 ellipsis               0.3.2     2021-04-29 [1] RSPM (R 4.1.0)
 evaluate               1.0.3     2025-01-10 [1] RSPM (R 4.1.0)
 farver                 2.1.2     2024-05-13 [1] RSPM (R 4.1.0)
 fastcluster          * 1.2.6     2024-01-12 [1] RSPM (R 4.1.0)
 fastmap                1.2.0     2024-05-15 [1] RSPM (R 4.1.0)
 foreach                1.5.2     2022-02-02 [1] RSPM (R 4.1.0)
 foreign                0.8-81    2020-12-22 [2] CRAN (R 4.1.2)
 Formula                1.2-5     2023-02-24 [1] RSPM (R 4.1.0)
 fs                     1.6.5     2024-10-30 [1] RSPM (R 4.1.0)
 genefilter             1.76.0    2021-10-26 [1] Bioconductor
 geneplotter            1.72.0    2021-10-26 [1] Bioconductor
 generics               0.1.3     2022-07-05 [1] RSPM (R 4.1.0)
 GenomeInfoDb         * 1.30.1    2022-01-30 [1] Bioconductor
 GenomeInfoDbData       1.2.7     2023-01-04 [1] Bioconductor
 GenomicRanges        * 1.46.1    2021-11-18 [1] Bioconductor
 getopt                 1.20.4    2023-10-01 [1] RSPM (R 4.1.0)
 ggbiplot             * 0.55      2023-01-22 [1] Github (vqv/ggbiplot@7325e88)
 ggplot2              * 3.5.1     2024-04-23 [1] RSPM (R 4.1.0)
 ggrepel              * 0.9.6     2024-09-07 [1] RSPM (R 4.1.0)
 ggside               * 0.3.1     2024-03-01 [1] RSPM (R 4.1.0)
 glue                   1.8.0     2024-09-30 [1] RSPM (R 4.1.0)
 GO.db                  3.14.0    2023-01-16 [1] Bioconductor
 gplots               * 3.2.0     2024-10-05 [1] RSPM (R 4.1.0)
 gridExtra              2.3       2017-09-09 [1] RSPM (R 4.1.0)
 gtable                 0.3.6     2024-10-25 [1] RSPM (R 4.1.0)
 gtools                 3.9.5     2023-11-20 [1] RSPM (R 4.1.0)
 Hmisc                  5.1-3     2024-05-28 [1] RSPM (R 4.1.0)
 htmlTable              2.4.3     2024-07-21 [1] RSPM (R 4.1.0)
 htmltools              0.5.8.1   2024-04-04 [1] RSPM (R 4.1.0)
 htmlwidgets            1.6.4     2023-12-06 [1] RSPM (R 4.1.0)
 httpuv                 1.6.15    2024-03-26 [1] RSPM (R 4.1.0)
 httr                   1.4.7     2023-08-15 [1] RSPM (R 4.1.0)
 impute                 1.68.0    2021-10-26 [1] Bioconductor
 IRanges              * 2.28.0    2021-10-26 [1] Bioconductor
 iterators              1.0.14    2022-02-05 [1] RSPM (R 4.1.0)
 KEGGREST               1.34.0    2021-10-26 [1] Bioconductor
 KernSmooth             2.23-20   2021-05-03 [2] CRAN (R 4.1.2)
 knitr                  1.50      2025-03-16 [1] RSPM (R 4.1.0)
 later                  1.4.1     2024-11-27 [1] RSPM (R 4.1.0)
 lattice              * 0.20-45   2021-09-22 [2] CRAN (R 4.1.2)
 lifecycle              1.0.4     2023-11-07 [1] RSPM (R 4.1.0)
 locfit                 1.5-9.12  2025-03-05 [1] RSPM (R 4.1.0)
 lpsymphony             1.22.0    2021-10-26 [1] Bioconductor (R 4.1.2)
 Maaslin2             * 1.8.0     2021-10-26 [1] Bioconductor
 magrittr               2.0.3     2022-03-30 [1] RSPM (R 4.1.0)
 MASS                   7.3-54    2021-05-03 [2] CRAN (R 4.1.2)
 Matrix                 1.5-3     2022-11-11 [1] RSPM (R 4.1.0)
 MatrixGenerics       * 1.6.0     2021-10-26 [1] Bioconductor
 matrixStats          * 1.5.0     2025-01-07 [1] RSPM (R 4.1.0)
 memoise                2.0.1     2021-11-26 [1] RSPM (R 4.1.0)
 mgcv                   1.8-38    2021-10-06 [2] CRAN (R 4.1.2)
 mime                   0.13      2025-03-17 [1] RSPM (R 4.1.0)
 miniUI                 0.1.1.1   2018-05-18 [1] RSPM (R 4.1.0)
 munsell                0.5.1     2024-04-01 [1] RSPM (R 4.1.0)
 mvtnorm                1.3-3     2025-01-10 [1] RSPM (R 4.1.0)
 nlme                   3.1-153   2021-09-07 [2] CRAN (R 4.1.2)
 nnet                   7.3-16    2021-05-03 [2] CRAN (R 4.1.2)
 optparse               1.7.5     2024-04-16 [1] RSPM (R 4.1.0)
 patchwork            * 1.3.0     2024-09-16 [1] RSPM (R 4.1.0)
 pcaPP                  2.0-5     2024-08-19 [1] RSPM (R 4.1.0)
 permute              * 0.9-7     2022-01-27 [1] RSPM (R 4.1.0)
 pillar                 1.10.1    2025-01-07 [1] RSPM (R 4.1.0)
 pkgbuild               1.4.6     2025-01-16 [1] RSPM (R 4.1.0)
 pkgconfig              2.0.3     2019-09-22 [1] RSPM (R 4.1.0)
 pkgload                1.4.0     2024-06-28 [1] RSPM (R 4.1.0)
 plyr                 * 1.8.9     2023-10-02 [1] RSPM (R 4.1.0)
 png                    0.1-8     2022-11-29 [1] RSPM
 preprocessCore         1.56.0    2021-10-26 [1] Bioconductor
 profvis                0.4.0     2024-09-20 [1] RSPM (R 4.1.0)
 promises               1.3.2     2024-11-28 [1] RSPM (R 4.1.0)
 purrr                  1.0.4     2025-02-05 [1] RSPM (R 4.1.0)
 R6                     2.6.1     2025-02-15 [1] RSPM (R 4.1.0)
 RColorBrewer           1.1-3     2022-04-03 [1] RSPM (R 4.1.0)
 Rcpp                   1.0.14    2025-01-12 [1] RSPM (R 4.1.0)
 RCurl                  1.98-1.17 2025-03-22 [1] RSPM (R 4.1.0)
 remotes                2.5.0     2024-03-17 [1] RSPM (R 4.1.0)
 reshape              * 0.8.9     2022-04-12 [1] RSPM (R 4.1.0)
 rlang                  1.1.5     2025-01-17 [1] RSPM (R 4.1.0)
 rmarkdown              2.29      2024-11-04 [1] RSPM (R 4.1.0)
 robustbase             0.99-4-1  2024-09-27 [1] RSPM (R 4.1.0)
 rpart                  4.1-15    2019-04-12 [2] CRAN (R 4.1.2)
 RSQLite                2.3.9     2024-12-03 [1] RSPM (R 4.1.0)
 rstudioapi             0.17.1    2024-10-22 [1] RSPM (R 4.1.0)
 S4Vectors            * 0.32.4    2022-03-24 [1] Bioconductor
 scales               * 1.3.0     2023-11-28 [1] RSPM (R 4.1.0)
 sessioninfo            1.2.3     2025-02-05 [1] RSPM (R 4.1.0)
 shiny                  1.10.0    2024-12-14 [1] RSPM (R 4.1.0)
 stringi                1.8.4     2024-05-06 [1] RSPM (R 4.1.0)
 stringr              * 1.5.1     2023-11-14 [1] RSPM (R 4.1.0)
 SummarizedExperiment * 1.24.0    2021-10-26 [1] Bioconductor
 survival               3.2-13    2021-08-24 [2] CRAN (R 4.1.2)
 tibble                 3.2.1     2023-03-20 [1] RSPM (R 4.1.0)
 tidyselect             1.2.1     2024-03-11 [1] RSPM (R 4.1.0)
 tximport             * 1.22.0    2021-10-26 [1] Bioconductor
 urlchecker             1.0.1     2021-11-30 [1] RSPM (R 4.1.0)
 usethis              * 3.1.0     2024-11-26 [1] RSPM (R 4.1.0)
 vctrs                  0.6.5     2023-12-01 [1] RSPM (R 4.1.0)
 vegan                * 2.6-10    2025-01-29 [1] RSPM (R 4.1.0)
 viridis              * 0.6.5     2024-01-29 [1] RSPM (R 4.1.0)
 viridisLite          * 0.4.2     2023-05-02 [1] RSPM (R 4.1.0)
 WGCNA                * 1.73      2024-09-18 [1] RSPM (R 4.1.2)
 withr                  3.0.2     2024-10-28 [1] RSPM (R 4.1.0)
 xfun                   0.51      2025-02-19 [1] RSPM (R 4.1.0)
 XML                    3.99-0.18 2025-01-01 [1] RSPM (R 4.1.0)
 xtable                 1.8-4     2019-04-21 [1] RSPM (R 4.1.0)
 XVector                0.34.0    2021-10-26 [1] Bioconductor
 zlibbioc               1.40.0    2021-10-26 [1] Bioconductor
```


