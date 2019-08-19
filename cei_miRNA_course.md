# Workshop

## Small non coding RNA

> 27 - 28 de septiembre en CICESE
>
> PhD Cei Abreu, LANGEBIO

#### Testing ...

Starts using virtual conection to the cluster via ssh

```bash
ssh rgomez@omica
# Password: ******
```

Set your work directory (Eg. $HOME/) and load the most current version of R available: `module load R-3.6.1` , Then start 

```R
ip <- installed.packages()[,c(1,3:4)]
rownames(ip) <- NULL
ip <- ip[is.na(ip[,3]),1:2,drop=FALSE]
print(ip, row.names=FALSE)
```

If local machine is used, please attend the follow instalation step:

**Installation**

```R
install.packages("BiocManager")

BiocManager::install()

.pkgs <- c("Rsamtools", "GenomicFeatures", "GenomicAlignments",
        "ggplot2", "pheatmap", "RColorBrewer","AnnotationDbi",
        "rtracklayer","GenomicRanges","edgeR","statmod","readr","biomaRt",
        "Biostrings","Biobase","knitr","affy","RCurl")

BiocManager::install(.pkgs)
```

Check Installation

```R

.inst <- .pkgs %in% installed.packages()

if(any(!.inst)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(.pkgs[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.pkgs), require, character.only = TRUE)
```

Session Info

```R
> sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server release 6.7 (Santiago)

Matrix products: default
BLAS:   /LUSTRE/apps/R-3.6.0/lib64/R/lib/libRblas.so
LAPACK: /LUSTRE/apps/R-3.6.0/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] RCurl_1.95-4.12             bitops_1.0-6
 [3] affy_1.62.0                 knitr_1.24
 [5] biomaRt_2.40.3              readr_1.3.1
 [7] statmod_1.4.32              edgeR_3.26.7
 [9] limma_3.40.6                rtracklayer_1.44.2
[11] RColorBrewer_1.1-2          pheatmap_1.0.12
[13] ggplot2_3.2.1               GenomicAlignments_1.20.1
[15] SummarizedExperiment_1.14.1 DelayedArray_0.10.0
[17] BiocParallel_1.18.1         matrixStats_0.54.0
[19] GenomicFeatures_1.36.4      AnnotationDbi_1.46.0
[21] Biobase_2.44.0              Rsamtools_2.0.0
[23] Biostrings_2.52.0           XVector_0.24.0
[25] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0
[27] IRanges_2.18.1              S4Vectors_0.22.0
[29] BiocGenerics_0.30.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2             locfit_1.5-9.1         lattice_0.20-38
 [4] prettyunits_1.0.2      assertthat_0.2.1       zeallot_0.1.0
 [7] digest_0.6.20          R6_2.4.0               backports_1.1.4
[10] RSQLite_2.1.2          httr_1.4.1             pillar_1.4.2
[13] zlibbioc_1.30.0        rlang_0.4.0            progress_1.2.2
[16] lazyeval_0.2.2         blob_1.2.0             Matrix_1.2-17
[19] preprocessCore_1.46.0  stringr_1.4.0          bit_1.1-14
[22] munsell_0.5.0          compiler_3.6.0         xfun_0.8
[25] pkgconfig_2.0.2        tibble_2.1.3           GenomeInfoDbData_1.2.1
[28] XML_3.98-1.20          crayon_1.3.4           withr_2.1.2
[31] grid_3.6.0             gtable_0.3.0           DBI_1.0.0
[34] magrittr_1.5           scales_1.0.0           stringi_1.4.3
[37] affyio_1.54.0          vctrs_0.2.0            tools_3.6.0
[40] bit64_0.9-7            hms_0.5.0              colorspace_1.4-1
[43] BiocManager_1.30.4     memoise_1.1.0
```



## **Linux tools** 

- samtools [http://www.htslib.org/](http://www.htslib.org/)
- bowtie1 [http://bowtie-bio.sourceforge.net/index.shtml](http://bowtie-bio.sourceforge.net/index.shtml)
- ShortStack https://github.com/MikeAxtell/ShortStack



```bash
# Para usarlo agregar esta trayectoria a la variable ambiental PATH:
export PATH=/LUSTRE/apps/bioinformatica/Shortstack:$PATH
Shortstack

#Este requiere de las otras aplicaciones que solicitaban, bowtie1 y samtools, del primero se encuentran las versiones 1.1.2 y 1.2.1.1 en las carpetas:
# <-- versión 1.1.2
export PATH=/LUSTRE/apps/bioinformatica/bowtie1:$PATH   

# <-- versión 1.2.1.1
export PATH=/LUSTRE/apps/bioinformatica/bowtie-1.2.1.1:$PATH  

# y samtools
export PATH=/LUSTRE/apps/bioinformatica/samtools-1.7/bin:$PATH

# Para usarlos en conjunto con Shortstack

export PATH=/LUSTRE/apps/bioinformatica/Shortstack:/LUSTRE/apps/bioinformatica/samtools-1.7/bin:/LUSTRE/apps/bioinformatica/bowtie-1.2.1.1:$PATH
```



Por instalar: `ViennaRNA`