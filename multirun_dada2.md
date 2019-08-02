Requerimos n formatos `RData`  y una lista de muestras en formato txt.

Este archivo debe contener el nombre de 1 muestra por linea. El formato del nombre es el usado en ciom_coi_dada2_multirun.R:``<crucero>_<estacion>_<tipo-muestra>`` 

Como generar el archivo para multi-run analisis en dada2:

1. Ir a directorio con los archivos fasta.gz
2. Extraer las muestras de un crucero (cambiar el primer 'ls' segun sea necesario)
3. Ajustar el nombre del archivo

```bash
 ls *fastq.gz | grep R1 | cut -f 1 -d '_' | cut -f 2,3,5 -d '-' > runX.samples
```

4. De ser necesario, abrir el archivo en un editor de texto, eliminar las muestras no desadas
5. Mover el runX.samples al directorio de analisis 

> **NOTA**: es necesario que el subjifo runX sea consistente con el subfijo del archivo RData; Por default, el sufijo de los objetos RData es *results.RData, debido a que estos se generan con el flujo de analisis del repositorio asv_coi, por tanto no es necesario incluir en la linea de comandos el nombre de los archivos debido a que estas extenciones (`result.RData` y `runX.samples` ) seran leidas por el script por si mismo.

6. Correr el siguiente script `cigom_coi_dada2_multirun.R`

```R
# cigom_coi_dada2_multirun.R
# Agosto 2019
#

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

if (is.na(args[1])) {
  threads <- TRUE
} else
  threads = as.numeric(args[1])

# ==============
## Checking and Load packages ----
# ==============
.cran_packages <- c("ggplot2", "GGally", "reshape2", "dplyr")
.bioc_packages <- c("dada2")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# +++++++++++++++

# ================================
## DADA2 Parameters in the config:
# ================================

# source("config.R")

date <- format(Sys.time(), "%Y%m%d")

out_prefix <- "multirun" # Output prefix (for plots and files)
run <- paste0("multirun_", date,"_COI")

# ================
# Outputs in `pwd`:
# ================
dir = getwd()
out_path <- file.path(dir, run)
system(command = paste0("mkdir -p ", out_path), intern = F)

## Merge seqtabs ----

blah <- list.files(pattern="results.RData", path=dir, full.names=TRUE)

foo <- list.files(pattern=".samples", path=dir, full.names=TRUE)

# Function to load RData and keep only one object
get_seqtab <- function(f, s) {
  # Obtener lista de muestras de interes
  mysamples <- read.table(s, header=F, stringsAsFactors = F)
  mysamples <- as.character(mysamples$V1)
  # Tal vez se puede aniadir sacar las columnas importantes (objeto 'out' y 'len_df') para la track.tsv y para histogram
  newenv <- new.env()
  load(file=f, env=newenv)
  f_seqtab <- newenv$seqtab
  # Subset 
  f_seqtab <- f_seqtab[which(rownames(f_seqtab) %in% mysamples),]
  f_seqtab <- f_seqtab[,colSums(f_seqtab) >= 1]
  rm(newenv)
  cat("\n", nrow(f_seqtab), "Samples and",ncol(f_seqtab), "number of ASVs selected to add in the output\n")
  return(f_seqtab)
}

# Load first seqtab
seqtab <- get_seqtab(blah[1], foo[1])
# Add other seqtabs
for(d in 2:length(blah)){
  d_seqtab <- get_seqtab(blah[d], foo[d])
  seqtab <- mergeSequenceTables(seqtab, d_seqtab)
}

# Report
cat("\nMerged seqtabs:\n")
dim(seqtab)

# Collapse no-mismatch
minOverlap <- 20

seqtab <- collapseNoMismatch(seqtab, minOverlap = minOverlap, orderBy = "abundance", verbose=TRUE)
# Report
cat("\nCollapsed no mismatch:\n")
dim(seqtab)


# Remove Quimeras ----
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = T)
cat("\nNo Bimera dimension:\n")
dim(seqtab.nochim)

# Check time
cat("\nTime after Bimera detection:\n")
system("date '+DATE: %Y-%m-%d%nTIME: %H:%M:%S'")

# Print Tables: ESVs
esv_out <- data.frame(seqtab.nochim)
esv_out$Sample <- rownames(seqtab.nochim)
esv_out <- esv_out[,c("Sample", colnames(seqtab.nochim))]
write.table(esv_out, file=paste0(out_path, "/", out_prefix, "_asv.tsv"), sep="\t", row.names = F, col.names = T)

# New distribution of Lengths
newlen_df <- data.frame(nchar(getSequences(seqtab.nochim)))
names(newlen_df) <- c("Length")
newlen_df$Process <- "Nochim"
cat("\nNew Length distribution:\n")
summary(newlen_df$Length)

# nice colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Prepare data and plot histogram
#lens_df <- rbind(len_df, newlen_df)
lens_df <- newlen_df
head(lens_df)
lens_plot <- ggplot(lens_df, aes(Length, color=Process)) + 
    geom_freqpoly(binwidth=1, size=1, alpha=0.7) +
    scale_color_manual(values = c("#999999", "#E7B800"))+ 
    labs(title=paste0(run, ". ESVs length distribution")) +
    theme_minimal() +
    theme(legend.position = c(0.05,0.95),
          legend.justification = c(0,1))

ggsave(paste0(out_prefix, "_lendistrib.jpeg"), plot=lens_plot, path=out_path)

# Proporcion de ESVs identificadas como quimeras
bimera_pct <- (1 - ncol(seqtab.nochim)/ncol(seqtab))*100
cat("\nASVs proportion identified as bimeras: ", bimera_pct, "%", "\n")

# Proporcion de reads que sobreviven
remain_reads <- (sum(seqtab.nochim)/sum(seqtab))*100
cat("\nSurvival ASVs: ", remain_reads, "%",  "\n")

# Track reads ----


merged <- rowSums(seqtab)
nochim <- rowSums(seqtab.nochim)
track <- data.frame(merged=merged, nochim=nochim)

# Print tables: Track
track_out <- data.frame(track)
track_out$Sample <- rownames(track)
track_out <- track_out[,c("Sample", colnames(track))]
write.table(track_out, file=paste0(out_path, "/", out_prefix, "_track.tsv"), sep="\t", row.names = F, col.names = T)

# Bar plot. Track
track_m <- melt(track_out, id.vars="Sample", variable.name = "Process", value.name = "Reads")
track_m$Process <- factor(track_m$Process)

# Color track
color_track <- function(track, s) {
  namef <- strsplit(s, "/")[[1]]
  namef <- namef[length(namef)]
  namef <- strsplit(namef, "_")[[1]][1]
  
  mysamples <- read.table(s, header=F, stringsAsFactors = F)
  mysamples <- as.character(mysamples$V1)
  
  #track[(rownames(track) %in% mysamples), "Factor"] <- namef
  track[(track$Sample %in% mysamples), "Factor"] <- namef
  return(track)
 }

for(d in 1:length(foo)){ track_m <- color_track(track_m, foo[d])}

# And plot
track_plot <- ggplot(track_m, aes(x=Sample, y=Reads, fill=Process)) +
    geom_col(position = position_identity(), width = 0.8, aes(fill = Process), alpha = 0.7) +
    scale_fill_brewer(palette = "BrBG", direction = 1) +
    labs(title=paste0(run, ". Reads count through processing"), x="") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 7), axis.text.y = element_text(size = 5)) + coord_flip()

track_plot <- track_plot + facet_grid(Factor ~ ., space = "free", scales = "free")

ggsave(paste0(out_prefix, "_trackplot.jpeg"), plot=track_plot, path=out_path)


# =================
# Save fasta file
# ================
# help in https://astrobiomike.github.io/amplicon/dada2_workflow_ex
asv_seqs <- colnames(seqtab.nochim)
n_asv <- dim(seqtab.nochim)[2]
asv_headers <- vector(n_asv, mode="character")

for (i in 1:n_asv) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))

write(asv_fasta,
            file=paste0(out_path, "/", 
                        out_prefix, "_ASVs.fasta")
                   )

# ============
# Save count table:
# ============ 
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)

write.table(asv_tab,
            file=paste0(out_path, "/", 
                        out_prefix, "_ASVs_count.table"), 
            sep="\t", 
            row.names = TRUE, 
            col.names = TRUE,
            quote=FALSE
                   )

# ============
# Save Image
# ============ 

#cat("\nSaving Outputs at\n", paste0(out_path, "/", out_path, ".results.RData"))
#save.image(file = paste0(out_path, "/", out_path, ".results.RData"))

cat("\nSaving Outputs at\n", paste0(out_path, "/multiresults.RData"))
save.image(file = paste0(out_path,"/multiresults.RData"))

# ======
# Finish
# ======

cat("\n**********\n**********\nI'm done here\nContinue w/ assignation step**********\n**********\n\n\n")

quit(save="no")
```



## Probando con crucero X05 de diferentes corridas

> **Ambientales**

```bash
ls 012*fastq.gz | grep R1 | cut -f 1 -d '_' | cut -f 2,3,5 -d '-'  > run012_X05.samples
#
ls 015*fastq.gz | grep R1 | cut -f 1 -d '_' | cut -f 2,3,5 -d '-'  > run015_X05.samples
```

Entonces

```bash
wc -l *samples
 # 6 run012_X05.samples
 # 27 run015_X05.samples
 # 33 total
```

Y corremos el siguiente script `sbatch.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=big_dada
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

#######################
# setting work_directory
# and exporting tools
#######################

module load R-3.5.0

echo 'Starting dada analysis within R!'
echo " "

Rscript --vanilla cigom_coi_dada2_multirun.R $SLURM_NPROCS
```

> **Ictioplancton**

```bash
ls 014*fastq.gz | grep R1 | cut -f 1 -d '_' | cut -f 2,3,5 -d '-'  > run014_X05.samples
#
ls 015*fastq.gz | grep R1 | cut -f 1 -d '_' | cut -f 2,3,5 -d '-'  > run015_X05.samples
# 
wc -l *samples

# 4 run014_X05.samples
# 39 run015_X05.samples
# 43 total

```

Entonces:

```bash
sbatch sbatch.sh
```

