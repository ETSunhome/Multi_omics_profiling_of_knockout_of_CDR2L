
# Multi-omics profiling reveals dysregulated ribosome biogenesis and impaired cell proliferation following knockout of CDR2L 

This repository holds all data and code to reproduce the main results, figures and supplementary data presented in:  
Solheim ET, Gerking Y, Kråkenes T, Herdlevær IV, Birkeland E, Totland C, Dick F, Vedeler CA (2024): “Multi-omics profiling reveals dysregulated ribosome biogenesis and impaired cell proliferation following knockout of CDR2L”. BMC Cancer. 2024;24(1):645. [doi: 10.1186/s12885-024-12399-z](10.1186/s12885-024-12399-z).

## Prerequisites

* R version 4.3.3.  
* All R packages used in the analysis are listed at the top of each Rmd file.  
* Several of the packages are bioconductor packages so it helps to have BiocManager installed:

```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(<"packagename">)
```

## How to..

### Re-run analysis from scratch

* The raw data consists of transcript-level counts in `./rawData/salmonOutput` and protein intensities in `./rawData/rawDataProteome.txt` and `./rawData/rawDataSecretome.txt`. See the Materials and Methods section of the published paper for details of how the datasets were created.  
* With these and the metadata in `./rawData/metadataTranscriptome`, `./rdsData/metadataProteome` and `./rdsData/metadataSecretome` all results can be reproduced.  
* `./preprocess.Rmd` Loads the raw data and performs differential expression analysis.  
* `./GSEA.Rmd` Performs gene set enrichment analysis (GSEA).   
* `./analysis.Rmd` Loads the results from the differential expression analysis and GSEA. Contains the code to reproduce all figures and tables presented in the paper.  

### Reproduce figures

* All figures created in R found in the paper can be reproduced with `./analysis.Rmd`.  
* The code loads `.rds` files needed for the analysis from `./rdsData/` and `./results`. These files are generated in `./preprocess.Rmd` and `./GSEA.Rmd`.  
* The figures are available in `./figs`.
* Some of the figures have been edited with graphics software afterwards and may therefore appear slightly different in the paper.


## Content

* `./preprocess.Rmd` Contains the code to preprocess data (transcript counts and protein intensities).  
* `./GSEA.Rmd` Contains the code to perform gene set enrichment analysis.  
* `./analysis.Rmd` Contains the code to reproduce figures and tables presented in the paper.  
* `./rawData/` Holds the raw data and metadata used in `./preprocess.Rmd`.   
    + `./rawData/rawDataProteome.txt` and `./rawData/rawDataSecretome.txt` Holds protein intensities for the proteome and secretome datasets, respectively.
    + `./rawData/metadataTranscriptome.txt`, `./rawData/metadataProteome.txt` and `./rawData/metadataSecretome.txt` Holds the sample metadata.  
    + `./rawData/salmonOut/` Holds output from Salmon for each sample in the analysis. Remember to unzip compressed files before running `./preprocess.Rmd`.  
* `./rdsData/` Holds .rds files created in `./preprocess.Rmd`.   
    + `./rdsData/countMatrixFiltered` Holds filtered, gene-level counts.  
    + `./rdsData/CPM` Holds filtered, gene-level counts in counts per million (CPM).  
    + `./rdsData/dds` Holds the DESeqDataSet (dds) object used for differential expression analysis.  
    + `./rdsData/proteins` Holds protein intensities for proteome dataset. Contaminants removed. 
    +  `./rdsData/proteins_filt` Holds protein intensities for proteome dataset. Filtered for missing values.  
    + `./rdsData/proteins_keep` Holds protein intensities for proteome dataset. Filtered for lowly expressed proteins. Missing values imputed using random forest.  
    + `./rdsData/m.RF` Holds matrix of protein intensities for proteome dataset. Missing values imputed using random forest.
    + `./rdsData/proteins_CCM` Holds protein intensities for secretome dataset. Contaminants removed.  
    +  `./rdsData/proteins_CCM_filt` Holds protein intensities for secretome dataset. Filtered for missing values.  
    + `./rdsData/proteins_CCM_keep` Holds protein intensities for secretome dataset. Filtered for lowly expressed proteins. Missing values imputed using random forest.  
    + `./rdsData/m_CCM.RF` Holds matrix of protein intensities for secretome dataset. Missing values imputed using random forest.  
* `./refData/` Holds reference data needed for the analysis.  
* `./data` Holds co-IP and proliferation data needed for the analysis.  
* `./figs` Holds the figures created in `./analysis.Rmd`.    
* `./results` Holds the results of the differential expression analysis and the gene set enrichment analysis.  
* `./cytoscape` Holds the files used as input for network analysis performed in Cytoscape. Files are created in `./analysis.Rmd`.  
* `./supplement` Holds the supplementary data. Several of the files are created in `./analysis.Rmd`.  


## Version info

```{r}
sessionInfo()
R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22631)

Matrix products: default

locale:
[1] C
system code page: 65001

time zone: Europe/Oslo
tzcode source: internal

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggpubr_0.6.0                plotrix_3.8-4               readxl_1.4.3                openxlsx_4.2.5.2            ComplexHeatmap_2.18.0       circlize_0.4.16            
 [7] WebGestaltR_0.4.6           UpSetR_1.4.0                cowplot_1.1.3               ggvenn_0.1.10               VennDiagram_1.7.3           futile.logger_1.4.3        
[13] reshape2_1.4.4              fgsea_1.28.0                lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
[19] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.0               tidyverse_2.0.0            
[25] edgeR_4.0.16                limma_3.58.1                missForest_1.5              DESeq2_1.42.0               SummarizedExperiment_1.32.0 Biobase_2.62.0             
[31] MatrixGenerics_1.14.0       GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1        
[37] tximport_1.30.0             matrixStats_1.2.0           data.table_1.15.2           biomaRt_2.58.2              rmarkdown_2.25             

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      jsonlite_1.8.8          rstudioapi_0.16.0       shape_1.4.6.1           magrittr_2.0.3          GlobalOptions_0.1.2     zlibbioc_1.48.0        
  [8] vctrs_0.6.5             memoise_2.0.1           RCurl_1.98-1.14         rstatix_0.7.2           htmltools_0.5.8.1       S4Arrays_1.2.0          progress_1.2.3         
 [15] itertools_0.1-3         lambda.r_1.2.4          curl_5.2.1              broom_1.0.5             cellranger_1.1.0        SparseArray_1.2.4       plyr_1.8.9             
 [22] futile.options_1.0.1    cachem_1.0.8            whisker_0.4.1           igraph_2.0.2            lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
 [29] Matrix_1.6-5            R6_2.5.1                fastmap_1.1.1           GenomeInfoDbData_1.2.11 clue_0.3-65             digest_0.6.34           colorspace_2.1-0       
 [36] AnnotationDbi_1.64.1    RSQLite_2.3.7           filelock_1.0.3          randomForest_4.7-1.1    timechange_0.3.0        fansi_1.0.6             httr_1.4.7             
 [43] abind_1.4-5             compiler_4.3.3          rngtools_1.5.2          bit64_4.0.5             withr_3.0.0             doParallel_1.0.17       backports_1.4.1        
 [50] BiocParallel_1.36.0     carData_3.0-5           DBI_1.2.2               apcluster_1.4.11        ggsignif_0.6.4          rappdirs_0.3.3          DelayedArray_0.28.0    
 [57] rjson_0.2.21            tools_4.3.3             zip_2.3.1               glue_1.7.0              cluster_2.1.6           generics_0.1.3          gtable_0.3.4           
 [64] tzdb_0.4.0              hms_1.1.3               xml2_1.3.6              car_3.1-2               utf8_1.2.4              XVector_0.42.0          ggrepel_0.9.5          
 [71] foreach_1.5.2           pillar_1.9.0            vroom_1.6.5             splines_4.3.3           BiocFileCache_2.10.2    lattice_0.22-5          bit_4.0.5              
 [78] tidyselect_1.2.1        locfit_1.5-9.9          Biostrings_2.70.3       knitr_1.45              gridExtra_2.3           svglite_2.1.3           xfun_0.42              
 [85] statmod_1.5.0           factoextra_1.0.7        stringi_1.8.3           yaml_2.3.8              evaluate_0.23           codetools_0.2-19        cli_3.6.2              
 [92] systemfonts_1.0.5       munsell_0.5.1           Rcpp_1.0.12             dbplyr_2.5.0            png_0.1-8               XML_3.99-0.16.1         parallel_4.3.3         
 [99] blob_1.2.4              prettyunits_1.2.0       doRNG_1.8.6             bitops_1.0-7            scales_1.3.0            crayon_1.5.2            GetoptLong_1.0.5       
[106] rlang_1.1.3             formatR_1.14            fastmatch_1.1-4         KEGGREST_1.42.0        
```