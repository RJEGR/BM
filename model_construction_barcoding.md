## Summary

**<u>Input file</u> **

- **(peces_bold.fasta):** <u>62,544</u>

**<u>Output file:</u>**

- **(peces_bold_0.2_cntrds_100_COI_bos_taurus_clean_vs_peces_bold.afa):** 62,562 -18  referencias
- **(mafft_clean_vs_peces_bold.afa):** 62,559 - 15 referencias
- **(coi_profile_vs_peces_bold_no_gaps.hmm.align):** 62529

Quick-run (Preparing model):

```bash
sbatch prep_profile.sh peces_bold.fasta 0.2 100

# Test the classification of the centroids
sbatch rdp_assign.sh peces_bold_0.2_cntrds_100_ab.fa 70 peces_bold.tax

# add sequence to the centroids and align
# 1
cat COI_bos_taurus.fasta peces_bold_0.2_cntrds_100_ab.fa > peces_bold_0.2_cntrds_100_COI_bos_taurus.fa
# 2
mafft --maxiterate 1000 --localpair  peces_bold_0.2_cntrds_100_COI_bos_taurus.fa > peces_bold_0.2_cntrds_100_COI_bos_taurus.afa

# Visualize aligment, remove misaligments and use as model to align multiple sequences with mafft or hmm

sbatch mafft.slurm peces_bold.fasta peces_bold_0.2_cntrds_100_COI_bos_taurus_clean.afa

```

# 1. Script description

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

