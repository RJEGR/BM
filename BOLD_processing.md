## Summary

**<u>Input file</u> **

- **(peces_bold.fasta):** <u>62,544</u>

**<u>Output file:</u>**

- **(peces_bold_0.2_cntrds_100_COI_bos_taurus_clean_vs_peces_bold.afa):** 62,562 -18  referencias
- **(mafft_clean_vs_peces_bold.afa):** 62,559 - 15 referencias
- **(coi_profile_vs_peces_bold_no_gaps.hmm.align):** 62529

# remover el nu

Directories:

**Mafft**

- /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/coi_alignment_fish/mafft

**HMM**

- /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/coi_alignment_fish/HMMER/rep_otu_bos_taurus

**Decipher**

- 

Preparing model by:

```bash
sbatch prep_profile.sh peces_bold.fasta 0.2 100

# Test the classification of the centroids
sbatch rdp_assign.sh peces_bold_0.2_cntrds_100_ab.fa 70 peces_bold.tax

# add sequence to the centroids and align
# 1
cat COI_bos_taurus.fasta peces_bold_0.2_cntrds_100_ab.fa > peces_bold_0.2_cntrds_100_COI_bos_taurus.fa
# 2
mafft --maxiterate 1000 --localpair  peces_bold_0.2_cntrds_100_COI_bos_taurus.fa > peces_bold_0.2_cntrds_100_COI_bos_taurus.afa

# Visualize aligment , remove misaligments and use as model to align multiple sequences with mafft or hmm

sbatch mafft.slurm peces_bold.fasta peces_bold_0.2_cntrds_100_COI_bos_taurus_clean.afa

```

# 1. Prepering a model

Copy and paste the follow script in a sbatch of name `prep_profile.sh`

```bash
#!/bin/bash
#
#SBATCH -J pprf
#SBATCH -p cicese
#SBATCH -t 6-00:00:00
#SBATCH -e %N.%j.err
#SBATCH -o %N.%j.out
#SBATCH -N 1
#SBATCH -n 12

# how to:
# sbacth prep_profile.sh seqs.fa 0.2 2

export PATH=/LUSTRE/apps/bioinformatica/hmmer-3.1b1/bin:$PATH

cd $SLURM_SUBMIT_DIR

rm *.tmp

seqs=$1
id=$2 # from 0 to 0.9, 0.2 (metaxa2 - build_db) or 0.3 (elbretch - primerMinner) recommended

min=$3 # from 1 to n 

usrt=${seqs%.*}_sorted.tmp
udrp=${usrt%.*}_derep.tmp
udrp_srt=${udrp%.*}_sorted.tmp
#cntrds=${seqs%.*}_${id}_cntrds.tmp
cntrds=${seqs%.*}_${id}_cnsns.tmp
acntrds=${cntrds%.*}_${min}_ab.fa

out_folder=tmp_files_${seqs%.*}

# 1. if gaps, remove it with:

#awk -f linearizeFasta.awk < $seqs | awk '{gsub(/[-, .]/, "", $3); print $1"\n"$3}' > ${seqs%.*}_ungapped.tmp

# or

clean_fasta.py -f $seqs


seqs=$(ls *_filtered.fasta)

# 2.
usearch  -sortbylength $seqs -output $usrt

# 3.
usearch -derep_fulllength $usrt -output $udrp -sizeout

usearch  -sortbylength $udrp -output $udrp_srt
# 4.
# usearch -cluster_smallmem $udrp_srt -id  $id -centroids $cntrds -usersort -sizeout
# -consout instead of -centroids
usearch -cluster_smallmem $udrp_srt -id  $id -consout $cntrds -usersort -sizeout

# save representative centroides
usearch -sortbysize $cntrds -minsize $min -output $acntrds

# OTUs distribution and size
grep '^>' $udrp | cut -d"=" -f2 | sed 's/;$//g'| sort | uniq -c | sort -k2,2 -n > derep_distr.txt
grep '^>' $cntrds | cut -d"=" -f2 | sed 's/;$//g'| sort | uniq -c | sort -k2,2 -n > centroids_distr.txt

# 5.
# High accuracy (for <~200 sequences x <~2,000 aa/nt):
mafft --maxiterate 1000 --localpair  $acntrds > ${acntrds%.*}.afa

mkdir -p $out_folder
mv *.tmp $out_folder

# visually inspect the alignments and remove misalignedsequences and gaps. Then construct profile and run hhmsearch. Ex:

# hmmbuild ${acntrds%.*}.hmm ${acntrds%.*}.afa

# hmmsearch --cpu $SLURM_NPROCS --tblout $tbl --domtblout $domtbl -A $out --acc ${acntrds%.*}.hmm $seqs

exit

```

b. If add additional sequence use:

```bash
# Ex
cat COI_bos_taurus.fasta peces_bold_0.2_cntrds_100_ab.fa > peces_bold_0.2_cntrds_100_COI_bos_taurus.fa

# High accuracy (for <~200 sequences x <~2,000 aa/nt):
mafft --maxiterate 1000 --localpair  peces_bold_0.2_cntrds_100_COI_bos_taurus.fa > peces_bold_0.2_cntrds_100_COI_bos_taurus.afa
```

# 2. Multiple sequences aligment using the model as reference

a. Run the multiple sequences aligment using the model as reference using `hmmsearch`:

```bash
#!/bin/bash
#
#SBATCH -J hmmer
#SBATCH -p cicese
#SBATCH -t 6-00:00:00
#SBATCH -e hmmer.%N.%j.err
#SBATCH -o hmmer.%N.%j.out
#SBATCH -N 2
#SBATCH -n 22

export PATH=/LUSTRE/apps/bioinformatica/hmmer-3.1b1/bin:$PATH
cd $SLURM_SUBMIT_DIR
#cd ./HMMER

prfl=$1
seqs=$2
tag=${prfl%.*}_vs_${seqs%.*}
out=${tag}.hmm
tbl=${tag}.tblout
domtbl=${tag}.domtblout

hmmsearch --cpu $SLURM_NPROCS --tblout $tbl --domtblout $domtbl -A $out --acc $prfl $seqs

./esl-reformat afa $out > ${out%.*}.align
```

b. Or Running the multiple sequence aligment using the model as reference using `mafft`, (reduce inter-genic gaps).  Save follow code in a script of name `mafft.slurm`.

```bash
#!/bin/bash
#
#SBATCH -J mafft
#SBATCH -p cicese
#SBATCH -t 6-00:00:00
#SBATCH -e %N.%j.err
#SBATCH -o %N.%j.out
#SBATCH -N 2
#SBATCH -n 22

ref=$2
seqs=$1
tag=${ref%.*}_vs_${seqs%.*}
out=${tag}.afa

cd $SLURM_SUBMIT_DIR

mafft --anysymbol --thread $SLURM_NPROCS --adjustdirection --keeplength --mapout --add $seqs --reorder $ref > $out

exit

```

check the aligment profile in R before doing the multiple aligment

```R
rm(list = ls())

# Test the model prepared with userch --> mafft

.bioc_packages <- c("Biostrings", "DECIPHER", "IRanges", "ggplot2")

# Load packages into session, and print package version
sapply(.bioc_packages, require, character.only = TRUE)

dir <- c('/Users/cigom/metagenomics/db/prep_model_ref/')

prpf <- dir(path = dir, pattern = 'afa', full.names = TRUE)

seqs <- readDNAStringSet(prpf, format="fasta")
BrowseSeqs(seqs)


# remove intergenic gapped sequences
# DMFH201-16;size=722;

seqs_adj <- AdjustAlignment(seqs)

cleanSeqs <- c('BAST1302-13;size=1620;', 'DMFH201-16;size=722;' )

seqs_adj <- seqs_adj[!names(seqs_adj) %in% cleanSeqs]

BrowseSeqs(seqs_adj)

seqs_adj <- subseq(seqs_adj, 609, width(seqs_adj))


outtag <- sub(".afa","",prpf[2])

writeXStringSet(seqs_adj, filepath = paste0(outtag, '_clean.afa'), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

# mandamos el alineamiento de secuencias 


```

Run the alignment of the model in R, **instead of mafft**

```R
rm(list = ls())

.bioc_packages <- c("Biostrings","DECIPHER") # "phangorn"

# Load packages into session, and print package version
sapply(.bioc_packages, require, character.only = TRUE)

# file <- dir(pattern="centroids_ab.fa")[1]
file <- dir(pattern="peces_bold_0.2_cnsns_100_COI_bos_taurus.fa")

dna <- readDNAStringSet(file, format="fasta")

# name "Vertebrate Mitochondrial" ,
SGC1 <- getGeneticCode("SGC1")

gT <- lapply(order(width(dna), decreasing=TRUE),
             function(x) {
               attr(x, "height") <- 0
               attr(x, "label") <- names(dna)[x]
               attr(x, "members") <- 1L
               attr(x, "leaf") <- TRUE
               x
             })


attr(gT, "height") <- 0.5
attr(gT, "members") <- length(dna)
class(gT) <- "dendrogram"
# use the guide tree as input for alignment


DNA <- AlignTranslation(dna,
                        guideTree=gT,
                        iterations=0,
                        refinements=0,
                        geneticCode = SGC1)

DNA <- AdjustAlignment(DNA)
BrowseSeqs(DNA)
# 

cleanSeqs <- c('centroid=BAST1302-13;seqs=1620;size=1620;')

DNA <- DNA[!names(DNA) %in% cleanSeqs]

writeXStringSet(DNA, filepath = paste0(file, ".decipher.afa"), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")


# CHECK CODON

aa <- translate(dna, genetic.code = SGC1)

BrowseSeqs(aa)

start <- subseq(aa, 1,1)

we <- width(aa)

ws <- we -5 

end <- subseq(aa, we, we)


```

Check the centroids classification vs the sequence classification

```bash
options(stringsAsFactors = FALSE)

library(tidyverse)

url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'
source(url)

centroid.file <- 'peces_bold_0.2_cntrds_100_ab_70_rdp.peces_bold_outdir/peces_bold_0.2_cntrds_100_ab.peces_bold.wang.taxonomy'

#centroid.file <- 'peces_bold_0.2_cnsns_100_COI_bos_taurus_70_rdp.peces_bold_outdir/peces_bold_0.2_cnsns_100_COI_bos_taurus.peces_bold.wang.tax.summary'

fish.file <- 'peces_bold.tax'

dim(centroid <- read_rdp(centroid.file))
dim(fish_bold <- read.delim(fish.file, row.names = 1, sep = ';', header = FALSE)) # 62544


# table(centroid[,2])
# table(fish_bold[,1])

length(rep_names <- names(table(centroid[,4]))) # 19
length(names <- names(table(fish_bold[,3]))) # 2712

dim(rep_fish_bold <- fish_bold[fish_bold[,3] %in% rep_names ,]) # 2807
dim(norep_fish_bold <- fish_bold[!fish_bold[,3] %in% rep_names,]) # 59737

# sanity check
identical(names(table(rep_fish_bold[,3])), rep_names) # TRUE
nrow(rep_fish_bold) + nrow(norep_fish_bold) == nrow(fish_bold) # TRUE

# 2807 lecturas subsisten la clasificacion de generos de los centroides

# 

library(data.table)

x <- data.table(z = 'Yes', table(rep_fish_bold[,3]))
y <- data.table(z = 'Not', table(fish_bold[,3]))

res <- rbind(x,y)
res <- res[order(N)]

# sanity check
length(rep_names) + length(names) == nrow(res) # TRUE

# TEST
# cumulative distribution of reads
par(mfrow=c(1,2))
plot(ecdf(res$N), main = "Cumulative distribution of genus", xlab='Number of sequences')

sample <- filter(res, N <= 500)

plot(ecdf(sample$N), main = "Cumulative distribution above the threshold", xlab='Number of sequences')

# make pct of abundance

pct <- function(x) {x / sum(x) * 100}

res[, pct := (pct(N))]

sum(res$pct) # 100

aggregate(res[,'pct'], by=list(res$z), sum)
# group data in fractions
res[(N == 1), group := "A"]
res[(N >= 2 & N <= 99), group := "B"]
res[(N >= 100 & N <= 499), group := "C"]
res[(N >= 500), group := "D"]

table(res$group)
length(table(res$group))

# 

library(ggalluvial)

ggplot(data = res,
       aes(axis1 = z, axis2 = group,
           y = pct)) +
  scale_x_discrete(limits = c("Centroid", "Group"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = z)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal(base_size=16) +
  labs(title = "", x = "Flow of data" , y = "Reads (%)")
  
  
quit(save = 'no')
```

Prepare taxa and fasta from data-input, filter coi_5p marker:

To detect diversity in the taxonomy

```R
options(stringsAsFactors = FALSE)
library(tidyverse)
tbl <- 'bovidae_bold_data.txt'
tbl <- read.delim(tbl)

# 
# check outliers ids if FALSE
length(unique(order(tbl[1]))) == nrow(tbl[1])
# Does all sequence have markercode? (TRUE)
sum(table(tbl$markercode)) == nrow(tbl)

# table(tbl$markercode)

marker <- c("COI-5P")

dim(tbl <- tbl[tbl$markercode %in% marker,])

seqs <- tbl$nucleotides
id <- tbl[1]

linage <- tbl[,c(10,12,14,16,18,20,22)]

# save fasta
id_ <- sub("^", ">", id$processid)
ws <- c(rbind(id_, seqs))

write(ws, file=paste0("bovidae.fa"))

# Only taxonomy

wt <- unite(linage, 1:ncol(linage), col='tax', sep = ";")
wt$tax <- sub("$", ";", wt$tax)
wt <- cbind(id, wt)

write.table(wt, file=paste0("bovidae.tax"), sep=" ", 
            row.names = FALSE, 
            col.names = FALSE,
            quote=FALSE)

# Plot a visualization

library(data.table)
taxtb <- data.table(table(linage$genus_name))
names(taxtb) <- c("rank", "n")
# cumulative distribution of reads
par(mfrow=c(1,2))
plot(ecdf(taxtb$n), main = "Cumulative distribution of genus", xlab='Number of sequences')
sample <- filter(taxtb, n <= 100)
plot(ecdf(sample$n), main = "Cumulative distribution above the threshold", xlab='Number of sequences')

#taxtb$n <- taxtb$n / sum(taxtb$n) *100

taxtb <- taxtb[order(-n), ]
taxtb[n <= 10, rank := "Others"]

library(ggpubr)

ggbarplot(taxtb, x = "rank", y = "n",
          x.text.angle = 90,
          ylab = "Number of sequence",
          xlab = "Genus",
          rotate = TRUE,
          ggtheme = theme_minimal()
          #facet.by = "Type"
          ) + theme(axis.text.y = element_text(hjust = 1, size = 7))



```

Las secuencias centroides encontrados en `otus_minsize.fa` pueden ser usados para construir el modelo de alineamiento de COI como se describe:

Per-alineamiento de centroides al modelo `COI_bos_taurus.fasta`

```bash

# con sbatch y srun no corre esta herramienta, averiguar y ejecutar!!!

hmmalign COI_bos_taurus.hmm peces_bold.subset.fasta

srun hmmscan --cpu 24 --domtblout TrinotatePFAM.out Pfam-A.hmm good.Trinity.fasta.transdecoder.pep

hmmscan --cpu $SLURM_NPROCS --domtblout domtblout.out $pfam $pep

```

And get the aligned sequence from file:

```bash
./esl-reformat afa coi_profile_vs_peces_bold_no_gaps.hmm > coi_profile_vs_peces_bold_no_gaps.hmm.align
```

The most important number here is the first one, the sequence E-value. The E-value is the expected number of false positives (nonhomologous sequences) that scored this well or better. The E-value is a measure of statistical significance. The lower the E-value, the more significant the hit. I typically consider sequences with E-values < 10−3 or so to be significant hits.

**So operationally:**

- if both E-values are significant (<< 1), the sequence is likely to be homologous to your query.

- if the full sequence E-value is significant but the single best domain E-value is not, the target sequence is probably a multidomain remote homolog; but be wary, and watch out for the case where
  it’s just a repetitive sequence.

- The fact that the two scores are slightly different is therefore telling you that there’s a small amount of probability (uncertainty) that the domain lies somewhat outside the envelope bounds that HMMER has selected

Scan large data-set (with errors)

```bash
#!/bin/bash
#SBATCH --job-name=boldBuild
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

#hmm=$1
#file=$2
#tag=${file%.}

wrap='hmmsearch --cpu 48 --tblout peces_bold_no_gaps.dblout --domtblout peces_bold_no_gaps.domtblout -A peces_bold_no_gaps.fa --acc mafft.hmm peces_bold_no_gaps.fasta'

echo $wrap | sh

exit

# test
hmmsearch --tblout dblout.out --domtblout domtblout.out -A peces_bold.subset.hmm --acc coi_profile.hmm peces_bold.subset.fasta

exit

sbatch -J hmmsearch -e hmmsearch.err -o hmmsearch.out -n 24 -N 1 -t 6-00:00 -mem=100GB \
-wrap="hmmsearch --cpu 24 dblout.out --domtblout domtblout.out -A peces_bold.subset.hmm --acc coi_profile.hmm peces_bold.subset.fasta"



```

asd

```bash
#nhmmer [options] <query hmmfile|alignfile> <target seqfile>

nhmmer COI_bos_taurus peces_bold.subset.fasta
```

## DNA Translation 

use vertebrate-mitocondrial

```bash

```



## Construct Folmer marker for BOLD data-base (not used)

#### Metaxa2 - conserve mode

...Metaxa2 has so far been strictly limited to operation on the SSU and LSU rRNA genes, preventing its use for other DNA barcodes. Yet, the capability of Metaxa2 to achieve high precision for its classifications while maintaining relatively high sensitivity would be highly desirable also for alternative barcoding markers, … presents an update to Metaxa2 itself, allowing the use of custom databases. We also introduce the Metaxa2 Database Builder (`metaxa2_dbb`)—a software tool that allows users to create customized databases from DNA sequences and their associated taxonomic affiliations.

In the `conserved mode`, on the other hand, the database builder will first extract the barcoding region from the input sequences using models built from a reference sequence provided and the Metaxa2 extractor. It will then align all the extracted sequences using MAFFT and determine the conservation of each position in the alignment. When the criteria for degree of conservation are met, all conserved regions are extracted individually and are then re-aligned separately using MAFFT. The re-aligned sequences are used to build hidden Markov models representing the conserved regions with HMMER. In this mode, the classification database will only consist of the extracted full-length sequences

In the `hybrid mode`, finally, the database builder will cluster the input sequences at 20% identity using USEARCH, and then proceed with the conserved mode approach on each cluster separately.

To install a database, first run the command “metaxa2_install_database” without any options. This will produce a list similar to this one:

```bash
# tuvimos problemas para establecer la conexion a las bases de datos debido a error en curl: error while loading shared libraries: libssl.so.1.0.0:
```



```bash
# remove gaps
awk -f linearizeFasta.awk <  peces_bold.fasta  | awk '{gsub(/[-, .]/, "", $3); print $1"\n"$3}'  > peces_bold_no_gaps.fasta

# get subset and remove gaps
awk -f linearizeFasta.awk <  peces_bold.fasta | head -n1000 | awk '{gsub(/[-, .]/, "", $3); print $1"\n"$3}'  > peces_bold.subset.fasta

# Get taxonomy
grep '^>' peces_bold.subset.fasta | sed 's/>//g' | column -t > peces_bold.subset.tax


# build metaxa database
metaxa2_dbb \
	-m peces_bold.subset.fasta \
	-t peces_bold.subset.tax \
	-g COI_bos_taurus.fasta \
	-d ./metaxa2_db
	
	--full_length 658
	--single_profile
	--mode conserved
	-d ./metaxa2_db
	
```

## BOLD Processing data base to rdp - mothur format

```bash
# ls data/
bold_Animals_Acanthocephala.trim      bold_Animals_Sipuncula.trim
bold_Animals_Annelida.trim            bold_Animals_Tardigrada.trim
bold_Animals_Arthropoda.trim          bold_Animals_Xenoturbellida.trim
bold_Animals_Brachiopoda.trim         bold_Fungi_Ascomycota.trim
bold_Animals_Bryozoa.trim             bold_Fungi_Basidiomycota.trim
bold_Animals_Chaetognatha.trim        bold_Fungi_Chytridiomycota.trim
bold_Animals_Chordata.trim            bold_Fungi_Glomeromycota.trim
bold_Animals_Cnidaria.trim            bold_Fungi_Myxomycota.trim
bold_Animals_Cycliophora.trim         bold_Fungi_Zygomycota.trim
bold_Animals_Echinodermata.trim       bold_Plants_Bryophyta.trim
bold_Animals_Gnathostomulida.trim     bold_Plants_Chlorophyta.trim
bold_Animals_Hemichordata.trim        bold_Plants_Lycopodiophyta.trim
bold_Animals_Mollusca.trim            bold_Plants_Magnoliophyta.trim
bold_Animals_Nematoda.trim            bold_Plants_Pinophyta.trim
bold_Animals_Nemertea.trim            bold_Plants_Pteridophyta.trim
bold_Animals_Onychophora.trim         bold_Plants_Rhodophyta.trim
bold_Animals_Platyhelminthes.trim     bold_Protists_Chlorarachniophyta.trim
bold_Animals_Porifera.trim            bold_Protists_Ciliophora.trim
bold_Animals_Priapulida.trim          bold_Protists_Heterokontophyta.trim
bold_Animals_Rotifera.trim            bold_Protists_Pyrrophycophyta.trim
```

Tenemos 50 bases con las siguientes caracteristica cada uno:

```BASH

# headers example (trimed file):

1. GBSP0794-06 # processid
2. 79455 # phylum_taxID
3. Cycliophora  # phylum_name   
4. 533072  # class_taxID
5. Eucycliophora  # class_name 
6. 533073  # order_taxID
7. Symbiida   # order_name     
8. 533074  # family_taxID
9. Symbiidae   # family_name    
10. NA   # subfamily_taxID        
11. --  # subfamily_name
12. 79456 #   genus_taxID
13. Symbion # genus_name
14. 79457  #  species_taxID
15. Symbion sp.  # species_name
16. COI-5P  # voucher_type
17. GGATCGTTACTAGGTGAC---GATC # nucleotides

```

Procesamos del siguiente modo:

```bash


# a) Tomamos el campo 1, que corresponde al processid, como etiqueta 

# b) Campo 3,5,7,9,13,15 de la Clasificacion taxonomica en el orden. ('Filo', 'Clase', 'Orden', 'Familia', 'Genero', 'Especie')

# c1) 
kingdom=*Protists
# c2) 
cut -f1,3,5,7,9,13,15 ${kingdom}trim | sed 's/ //g' > taxonomy.${kingdom}.csv 
# replace ' ' for ''
Rscript --vanilla bold_public_process_for_RDP.R taxonomy.${kingdom}.csv FALSE




```

El siguiente script identifica el primer campo (columna, como el identificador de la especie [processid]); así como la última posición (ie. última columna) a la secuencia de origen [voucher_type];Si existen identificadores repetidos, el script re-etiquetara con el sufijo equivalente a las n repeticiones del identificador separado por un guión bajo “_”. El resto de las columnas son implementadas para formatear la taxonomía para el algoritmo de mothur (classify.seqs). El nombre del archivo (Ej. taxonomy.Animals.csv) es utilizado para etiquetar el segundo nivel taxonómico (Ej. root;rank1, donde rank1 pertenece al nivel taxonómico reino) por lo que es necesario nombrar el archivo de entrada con la etiqueta correspondiente. Los campos de la asignación taxonómica vacíos son rellenados con valores NA, y finalmente se retienen sólo los taxones con resolución hasta especie. 



El script puede ser descargado de [bold_public_process_for_RDP.R](https://raw.githubusercontent.com/RJEGR/metagenomics/master/bold_public_process_for_RDP.R)

> Script:

```r
#!/usr/bin/env Rscript

# El siguiente script identifica el primer campo (columna, como el identificador de la especie [processid]); 
# así como la última posición (ie. última columna) a la secuencia de origen [voucher_type];
# Si existen identificadores repetidos, el script re-etiquetara con el sufijo equivalente a las n repeticiones del identificador separado por un guión bajo “_”. 
# El resto de las columnas son implementadas para formatear la taxonomía para el algoritmo de mothur (classify.seqs). 
# El nombre del archivo (Ej. taxonomy.Animals.csv) es utilizado para etiquetar el segundo nivel taxonómico (Ej. root;rank1, donde rank1 pertenece al nivel taxonómico reino) por lo que es necesario nombrar el archivo de entrada con la etiqueta correspondiente. 
# Los campos de la asignación taxonómica vacíos son rellenados con valores NA, y finalmente se retienen sólo los taxones con resolución hasta especie. 

# How tu run:
# bold_Plants_Chlorophyta.trim, bold_Plants_Magnoliophyta.trim, bold_Plants_Bryophyta.trim
# ~$ kingdom='Plants'
# ~$ cut -f1,3,5,7,9,13,15,17 *${kingdom}*trim | sed 's/ /_/g' > taxonomy.${kingdom}.csv 
# Rscript --vanilla bold_public_process_for_RDP.R taxonomy.${kingdom}.csv FALSE

args = commandArgs(trailingOnly=TRUE)

cat('1. Reading file at time:', (d <- format(Sys.time(), "%H:%M:%S")), '\n')


file = args[1]

# By default keep species is TRUE
if (is.na(args[2])) {
    keep_spp <- TRUE
} else
keep_spp = args[2]

kingdom <- strsplit(file, '[.]')[[1]][2]
root.rank <- paste0('root;Eukaryota',';', kingdom)
colnames <- c('Id', 'Filo', 'Clase', 'Orden', 'Familia', 'Genero', 'Especie', 'Secuencia')
# read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
x <- read.csv(file, header=FALSE, sep='\t', na.strings=c("","NA"), stringsAsFactors=FALSE)
names(x) <- colnames

cat('2. Keeping only the Specie resolution ...\n')


if (keep_spp) { x <- x[complete.cases(x),] }


Id <- make.unique(as.vector(x[,1]), sep = "_")

cat('3. Adding root and ...\n')

y <- cbind(root.rank, (x[-c(1,ncol(x))]))


cat('4. Formating taxonomy ...\n')

library(tidyr)

save <- cbind(Id, unite(y, sep = ";", remove = TRUE, col = 'Taxonomy'))

save$Taxonomy <- sapply(save$Taxonomy, 
                        function(x){gsub(pattern = "$",
                        replacement = ";", x)})


cat('5. Formating fasta ...\n')

n_seqs <- length(Id)
seq_headers <- vector(n_seqs, mode="character")

for (i in 1:n_seqs) {
  seq_headers[i] <- paste(">", Id[i], sep = "")
}

fasta <- c(rbind(seq_headers, x[,ncol(x)]))

cat('6. Writing outputs ... \n')
# 1
write.table(save, file = paste0("Bold.",kingdom, ".tax"), append = FALSE, quote = FALSE, sep = " ",
            na = "NA", row.names = FALSE,
            col.names = FALSE)
# 2
write(fasta,
            file=paste0("Bold.",
                          kingdom, 
                          ".fasta"))

cat('7. DONE!\n')

difftime(strptime(format(Sys.time(), "%H:%M:%S"), format = "%H:%M:%S"), 
         strptime(d, format = "%H:%M:%S"),units="auto")

quit(save ='no')


```

Concatenamos 

```BASH
cat ./*/Bold.*.fasta > BOLD_public.ALL.fasta
cat ./*/Bold.*.tax > BOLD_public.ALL.tax
# or
cat ./*/Bold.*.tax > BOLD_public_species.tax
cat ./*/Bold.*.fasta > BOLD_public_species.fasta
```



Upload to cluster:

```bash
scp BOLD_public.ALL* rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/BOLD
```

and run a test:

```bash
sbatch rdp_assign.sh run012_relax_ASVs.fasta 99 BOLD_public.ALL.tax
```

> rdp_assign.sh script

```bash
#!/bin/bash
#SBATCH --job-name=rdp
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

# Variables

fasta=$1 # sequence to assign
boots=$2 # the interger of bootstrap to use in classify.seqs
DB_TAX=$3 # may be the taxa file, by default the script will be search for homolgy fasta name
DB_REF=${DB_TAX%.*}.fasta
outdir=${fasta%.*}_${boots}_rdp.${DB_TAX%.*}_outdir

## database
DB=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/BOLD/
#DB_REF="BOLD_public_species_v02.fasta"
#DB_TAX="BOLD_public_species_v02.tax"

mothur="/LUSTRE/bioinformatica_data/genomica_funcional/bin/mothur/mothur"

echo "Asignacion taxonomica de: $fasta"
echo "DB: $DB_REF and $DB_TAX"
echo "bootstrap cutoff: $boots"

cd $SLURM_SUBMIT_DIR

echo "Starting ... !"

$mothur "#system(mkdir -p $outdir);set.dir(output=$outdir, tempdefault=$DB);summary.seqs(fasta=$fasta, processors=$SLURM_NPROCS);classify.seqs(fasta=current, reference=$DB_REF, taxonomy=$DB_TAX, iters=1000, cutoff=$boots);get.current();quit()"

exit
```



test mothur single end

```bash
merge.files(input=cigom.files, output=cigom.fasta)

fastq.info(file=cigom.file)
fastq.info(fastq=016-X07-A3-18S-AMB_S1_L001_R1_001.fastq.gz)

```
