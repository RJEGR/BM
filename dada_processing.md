# Introduction to dada2

>  Benjamin J Callahan, Joey McMurdie, Susan Holmes

The investigation of environmental microbial communities and microbiomes has been revolutionized by the development of high-throughput amplicon sequencing. In amplicon sequencing a particular genetic locus, for example the 16S rRNA gene in bacteria, is amplified from DNA extracted from the community of interest, and then sequenced on a next-generation sequencing platform. This technique removes the need to culture microbes in order to detect their presence, and cost-effectively provides a deep census of a microbial community.

However, the process of amplicon sequencing introduces errors into the sequencing data, and these errors severely complicate the interpretation of the results. DADA2 implements a novel algorithm that models the errors introduced during amplicon sequencing, and uses that error model to infer the true sample composition. DADA2 replaces the traditional “OTU-picking” step in amplicon sequencing workflows, producing instead higher-resolution tables of amplicon sequence variants (ASVs).

As seen in [the paper introducing DADA2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/) and in [further benchmarking](http://benjjneb.github.io/dada2/R/SotA.html), the DADA2 method is more sensitive and specific than traditional OTU methods: DADA2 detects real biological variation missed by OTU methods while outputting fewer spurious sequences. [Another recent paper](https://www.nature.com/ismej/journal/vaop/ncurrent/full/ismej2017119a.html)describes how replacing OTUs with ASVs improves the precision, comprehensiveness and reproducibility of marker-gene data analysys.

## Preparing dataset

Lets first sampling a small subset

```bash
mkdir SAMPLING
for i in $(ls *AMB*R1* | shuf -n7); do echo $i | xargs -I {} cp {} SAMPLING; done
# and R2
cd SAMPLING
for i in $(ls *AMB*R1*); do echo ${i/R1/R2} | xargs -I {} cp ../{} .; done
```

> wd=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/data/run09_20180511_COI/SAMPLING

The libraries to use are:

```bash
50M 09-X04-A2-COI-AMB_S1_L001_R1_001.fastq.gz
56M 09-X04-A2-COI-AMB_S1_L001_R2_001.fastq.gz
40M 09-X04-A8-COI-AMB_S6_L001_R1_001.fastq.gz
46M 09-X04-A8-COI-AMB_S6_L001_R2_001.fastq.gz
46M 09-X04-C25-COI-AMB_S19_L001_R1_001.fastq.gz
50M 09-X04-C25-COI-AMB_S19_L001_R2_001.fastq.gz
42M 09-X04-D29F-COI-AMB_S23_L001_R1_001.fastq.gz
49M 09-X04-D29F-COI-AMB_S23_L001_R2_001.fastq.gz
30M 09-X04-E34-COI-AMB_S27_L001_R1_001.fastq.gz
34M 09-X04-E34-COI-AMB_S27_L001_R2_001.fastq.gz
43M 09-X04-G41-COI-AMB_S33_L001_R1_001.fastq.gz
47M 09-X04-G41-COI-AMB_S33_L001_R2_001.fastq.gz
44M 09-X04-Y1A-COI-AMB_S40_L001_R1_001.fastq.gz
50M 09-X04-Y1A-COI-AMB_S40_L001_R2_001.fastq.gz
```



## **Load data and define parameters**

Cargamos el modulo donde tenemos los datos`module load R-3.5.0`

1) Ruta de las bases de datos

```R

coi_taxa <- "/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/midori_long_coi_taxa.fa.gz" #omica

coi_spp <- "/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/midori_uniq_coi_taxa_spp.fa.gz" #omica

# Directorio de trabajo
dir <- getwd()
coi_path <- dir

```

2) Y parametros importantes

```R
# Do Filtering
do_derep <- TRUE
# Do Assign Taxonomy
do_at <- TRUE

# Run
run <- "run09_20180511_COI"

# Input prefix (to be trimmed)
in_prefix <- "09-X04-"

# Output prefix (for plots and files)
n_test <- "test04"

out_prefix <- paste0("run09_", n_test)
# Sample type, one of the following
sample_type <- "" # Todas las muestras
#sample_type <- "AMB.*" # Solo muestras ambientales
#sample_type <- "ICT.*" # Solo muestras de Ictio

## DADA2 Parameters: 
# FilterAndTrim (see ?filterAndTrim)
fat_trunclen <- c(230,170)
fat_trunQ <- 15
fat_maxee <- c(2,2)
fat_minlen <- 120
# Assign Taxonomy minimum bootstrap (see ?assignTaxonomy)
at_minboot <- 80

# Oficial Levels
TL <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

```

3) Cargamos paquetes o instalamos:

```R
# ==============
## Load packages ----
# ==============
if (!require('dada2')) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
        BiocManager::install("dada2", version = "3.8")
}  else
    if (!require("ggplot2")) { 
        install.packages('ggplot2', dep=TRUE, repos='http://cran.us.r-project.org') 
    } else
    if (!require("dplyr")) {
        install.packages('dplyr', dep=TRUE, repos='http://cran.us.r-project.org') 
    } else
    if (!require("reshape2")) {
        install.packages('reshape2', dep=TRUE, repos='http://cran.us.r-project.org')
     }
```

4) ...

```R

## DADA2 pipeline ----

# Directorio con archivos
system("mkdir data")
system("mv *gz ./data")
data_path <- paste0(coi_path, "/data/")     #omica
cat("\nFiles available in input dir:", data_path, "\n")
list.files(data_path)

# Directorio para outputs (plots and tables)
out_path <- file.path(coi_path, "dada2_esv", out_prefix) #omica
system(command = paste0("mkdir -p ", out_path), intern = F)

# Crear nombres de archivos
fnFs <- sort(list.files(data_path, pattern = paste0(sample_type, "_R1_001.fastq.gz"), full.names = T))
cat("\nForward read files to process:\n", fnFs, "\n")
fnRs <- sort(list.files(data_path, pattern = paste0(sample_type, "_R2_001.fastq.gz"), full.names = T))
cat("\nReverse read files to process:\n", fnRs, "\n")

sample.names <- sapply(strsplit(basename(fnFs), "_"), '[', 1)
sample.names <- sub(in_prefix, "", sample.names) #
cat("\nSample names:\n", sample.names, "\n")

```

## Filtering and trimming

```R
fnFsQC <- plotQualityProfile(fnFs[1:4])
ggsave(paste0(out_prefix, "_FsQC.jpeg"), plot=fnFsQC, path=out_path, width = 10, height = 7, units = "in")
rm(fnFsQC)
fnRsQC <- plotQualityProfile(fnRs[1:4])
ggsave(paste0(out_prefix, "_RsQC.jpeg"), plot=fnRsQC, path=out_path, width = 10, height = 7, units = "in")
rm(fnRsQC)
```

Ahora visualizamos el perfil de calidad de las lecturas en algunas muestras:![](/Users/cigom/metagenomics/COI/SAMPLING/dada2_esv/run09_test04/run09_test04_RsQC.jpeg)

En escala de colores grises (mapa de calor) se muestra la **frecuencia** del puntuaje de calidad para cada base. La **media** del puntuaje de calidad en cada posicion se muestra en la linea verde, mientras que los **cuartiles** de la distribucion del puntuaje de calidad en lineas naranjas. La línea roja nos indica la proporción escalada de lecturas que se extienden hasta al menos esa posición.

#### Lecturas Forward

En general las lecturas forward son de buena calidad, recomendamos recortar los últimos nucleótidos para evitar errores menos controlados que puedan surgir allí. Estos perfiles de calidad no sugieren que se necesite ningún recorte adicional. Truncaremos las lecturas hacia delante en la posición 230 (recortando los últimos 80 bases nucleótidos).

![](/Users/cigom/metagenomics/COI/SAMPLING/dada2_esv/run09_test04/run09_test04_FsQC.jpeg)

#### Lecturas Reverse

Las lecturas *reverse* suelen ser de una calidad significativamente mas baja, especialmente al final, lo cual es común en la secuenciación de Illumina. DADA2 incorpora información de calidad en su modelo de error, lo que hace que el algoritmo sea robusto a una secuencia de calidad más baja, pero el recorte a medida que se desplome la calidad promedio mejorará la sensibilidad del algoritmo a variantes de secuencia raras. Sobre la base de estos perfiles, truncaremos las lecturas inversas en la posición 160 donde la distribución de calidad falla.

The reverse reads are of significantly worse quality, especially at the end, which is common in Illumina sequencing. This isn’t too worrisome, as DADA2 incorporates quality information into its error model which makes the algorithm [robust to lower quality sequence](https://twitter.com/bejcal/status/771010634074820608), but trimming as the average qualities crash will improve the algorithm’s sensitivity to rare sequence variants. Based on these profiles, we will truncate the reverse reads at position 160 where the quality distribution crashes.

```R
Filter and Trim Parameters:
        truncLen:  230 170  Reads shorter than this bases are discarded
        truncQ:  15  Quality score less than or equal
        maxEE:  2 2  expected errors will be discarded
        minLen:  120  Remove reads with length less than minLen
        rm.phix:  TRUE  discard reads that match against the phiX genome

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = fat_trunclen, 
                     truncQ = fat_trunQ, 
                     maxEE= fat_maxee, 
                     minLen = fat_minlen, 
                     trimLeft = 26, maxN = 0, 
                     rm.phix = rm_phix, compress = T, multithread = T)


```

> We’ll use standard filtering parameters: `maxN=0` (DADA2 requires no Ns), `truncQ=2`, `rm.phix=TRUE` and `maxEE=2`. The `maxEE`parameter sets the maximum number of “expected errors” allowed in a read, which is [a better filter than simply averaging quality scores](https://academic.oup.com/bioinformatics/article/31/21/3476/194979).

Y como resultado:

![](/Users/cigom/metagenomics/COI/SAMPLING/dada2_esv/run09_test04/run09_test04_FsQC_filt.jpeg)

## Dada algorithm

The DADA2 algorithm makes use of a parametric error model (`err`) and every amplicon dataset has a different set of error rates. The `learnErrors` method learns this error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

 ```R
errF <- learnErrors(filtFs, MAX_CONSIST = 20, multithread = T)
errR <- learnErrors(filtRs, MAX_CONSIST = 20, multithread = T)
 ```



Outlines:

1. **Load data and parameters**

2. **Filtering and trimming**

   a. Plot quality and trimming

   b. Trimming

   c. Tracking data

3. **Dada algorithm** `threads selection`

   a. Error (learn Error)

   b. De-replicate

   c. dada algorithm 

   d. Merge Pairs

4. **Taxonomy assignation**

   a. At_min_boot `$variable config`

   b. Print taxa

   c. make.biom (or Rdata)

   d. add tree (usearch or phantom)

5. **Phyloseq.Rmd**

   a. Load data and parameters

   b. Plot quality and trimming

   c. Trimming

   d. Tracking data

# Run it from the cluster

Descargar flujo desde

`git clone https://github.com/RJEGR/dada2.git`

El script desde el cluster se manda dela siguiente manera, 

`sbatch dada2_coi.sh TESTNAME DATASET_PATH`

Ejemplo de resultados finales estaran aqui:

`/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/SCRIPTING/TEST_FULL_DATASET`

