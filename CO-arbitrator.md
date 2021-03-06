# CO-arbitrator

> Heller P. et al 2018

Procesar el formato de base de datos que los autores de CO-arbitrator desarrollaron:

Fasta file containing COI nucleotides records retrieved from GenBank by CO-ARBitrator. Deflines are delimited by double underscores, and consist of an accession number, followed by a binomial identification, followed by taxonomy. Values within the binomial and the taxonomy are delimited by semicolons.

_The binomial identification and the taxonomy are internally delimited by semicolons._

> La nomenclatura es la del GenBank

```bash
>NC_004385__Melanotaenia;lacustris__Eukaryota;Metazoa;Chordata;Craniata;Vertebrata;Euteleostomi;Actinopterygii;Neopterygii;Teleostei;Neoteleostei;Acanthomorphata;Ovalentaria;Atherinomorphae;Atheriniformes;Melanotaeniidae;Melanotaenia.
GTGATAATTACACGTTGATTCTTCTCTACTAATCACAAAGACATT ...
#
KJ940174__Undinula;vulgaris__Eukaryota;Arthropoda;Maxillopoda;Calanoida;Calanidae;Undinula;Undinulavulgaris #<-- el nivel especie esta unido
ATAATTGGAACAGGTTTAAGAATAATTATTCGATTAGAATTAGGACAAGCTGGATCATTAATTGGAGATGATCAGATTTATAATGTAGT...
```

El archivo `Coarbitrator_COI_nuc.fa` contiene 1,043,596 secuencias:

Procesamos:

```bash
# Obtenemos el archivo fasta, unicamente con el numero de acceso
file=Coarbitrator_COI_nuc.fa

awk '/^>/{gsub(/__/, " "); print ""$1; next}{print}' $file > ${file%.*}.fasta

# Ordenamos el identificador binomial (hace referencia a la especie) y el linaje taxonomico en la posicion rank K,P,C,O,F,G,S

# 1 BINOMIAL IDS
awk '/^>/{gsub(/__/, " "); print $1"\t"$2;}' ${file} | sed 's/;/ /g' | sed 's/^>//g' > ${file%.*}.binomial.ids

# 2 TAXA FILE
awk '/^>/{gsub(/__/, " "); print $1"\t"$3}' $file | sed 's/^>//g' | sed 's/.$/;/g' > ${file%.*}.tax
# extra.
# Count the number of max levels (between 2-19) - full lineage / abbreviate lineage
awk -F';' '{print (NF ? NF+0 : 0)}' Coarbitrator_COI_nuc.tax | sort | uniq -c | sort -n -k2,2
# al articulo reporta que 461 secuencias tuvieron asignacion nivel-genero y 199,119 identificaciones de mayor nivel
 112 3
 284 4
2452 5
54100 6
490819 7
5463 8
7072 9
18497 10 <-- sigue siendo genero
19574 11
17374 12
48731 13
121591 14
102675 15
116609 16
32780 17
4871 18
 564 19
  28 20
```

El archivo `Coarbitrator_COI_aa.faa` contiene 1,043,739,  (en el articulo determinan 1,054,973 secuencias proteicas [>= 95 aa] del metazoan-genBank clasificadas satisfactoriamente )

procesamos:

```bash
# Obtenemos el archivo fasta, unicamente con el numero de acceso
file=Coarbitrator_COI_aa.faa
awk '/^>/{gsub(/__/, " "); print ""$1; next}{print}' $file > ${file%.*}.fasta

# Ordenamos el identificador binomial (hace referencia a la especie) y el linaje taxonomico en la posicion rank K,P,C,O,F,G,S

# 1 BINOMIAL IDs
awk '/^>/{gsub(/__/, " "); print $1"\t"$2;}' ${file} | sed 's/;/ /g' | sed 's/^>//g' > ${file%.*}.binomial.ids

# 2 TAXA FILE
awk '/^>/{gsub(/__/, " "); print $1"\t"$3}' $file | sed 's/^>//g' | sed 's/.$/;/g' > ${file%.*}.tax


# extra.

# Count the number of max levels (between 2-19) - full lineage / abbreviate lineage
awk -F';' '{print (NF ? NF+0 : 0)}' ${file%.*}.tax | sort | uniq -c | sort -n -k2,2
# al articulo reporta que 461 secuencias tuvieron asignacion nivel-genero y 199,119 identificaciones de mayor nivel

 112 3
 291 4
2452 5
54101 6
490822 7
5482 8
7103 9
18504 10
19586 11
17386 12
48733 13
121606 14
102683 15
116622 16
32791 17
4872 18
 565 19
  28 20
```

Entonces, subimos al cluster y probamos base de datos



```bash
scp -r /Users/cigom/metagenomics/db/co-arbitrator rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs
```



### rdp classifier

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

fasta=$1
boots=$2
outdir=${fasta%.*}_${boots}_rdp_outdir

## CO-arbitrator
DB=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/co-arbitrator
DB_REF="Coarbitrator_COI_nuc.fasta"
DB_TAX="Coarbitrator_COI_nuc.tax"

mothur="/LUSTRE/bioinformatica_data/genomica_funcional/bin/mothur/mothur"

echo "Asignacion taxonomica de: $fasta"
echo "DB: $DB_REF and $DB_TAX"
echo "bootstrap cutoff: $boots"

cd $SLURM_SUBMIT_DIR

echo "Starting ... !"

$mothur "#system(mkdir -p $outdir);set.dir(output=$outdir, tempdefault=$DB);summary.seqs(fasta=$fasta, processors=$SLURM_NPROCS);classify.seqs(fasta=current, reference=$DB_REF, taxonomy=$DB_TAX, iters=1000, cutoff=$boots);get.current();quit()"

exit

```

Entonces

```bash
sbatch rdp_assign.sh run014_t2_ASVs.fasta 99 Coarbitrator_COI_nuc_curated.tax
```

Algunos **errores** en la base de datos debido a: 

- Secuencias no anadidas (1011) o rellenadas el campo con el termino 'none'

```bash
Generating search database...    [WARNING]: We found more than 25% of the bases in sequence complement(NC_008833 to be ambiguous...
#
grep  -A1'NC_008833' Coarbitrator_COI_nuc.fa

>complement(NC_008833__Placozoan;sp.;BZ49__Eukaryota;Metazoa;Placozoa;unclassifiedPlacozoa.
none
# 
grep -c 'none' Coarbitrator_COI_nuc.fa
1011
```

a. Identificadores repetidos.

b. Identificadores sin identificacion ('NONE')

```bash
awk '{print $1}' Coarbitrator_COI_nuc.tax | sort | uniq -c | sort -n -k2,2 | tail
  # a ---> 2 NC_006354
  # 2 NC_016423
  # 2 NC_016463
  # 2 NONE <---- b
  # 3 HE964902
```

Trabajaremos con R para limpiar los errores de los puntos a y b:

```R
#!/usr/bin/env Rscript
# reviewed version from ~/Documents/GitHub/metagenomics/bold_public_process_for_RDP.R

.bioc_packages <- c("Biostrings", "IRanges", "tidyr", "zoo")

# 2.
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(.bioc_packages[!.inst], ask = F)
}

sapply(c(.bioc_packages), require, character.only = TRUE)

#

path_db <- '/Users/cigom/metagenomics/db/co-arbitrator'
setwd(path_db)

fasta.file <- "Coarbitrator_COI_nuc.fasta"
tax.file <- "Coarbitrator_COI_nuc.tax"
binomail.file <- "Coarbitrator_COI_nuc.binomial.ids"

fasta.obj <- readDNAStringSet(fasta.file)
seqs <- as.data.frame(fasta.obj)$x

#   reading FASTA file Coarbitrator_COI_nuc.fasta: ignored 56 invalid one-letter sequence codes

taxa.obj <- read.csv(tax.file, header=FALSE, sep='\t', na.strings=c("","NA"), stringsAsFactors = FALSE)
# 
bim.obj <- read.csv(binomail.file, header=FALSE, sep='\t', stringsAsFactors = FALSE)

# anadimos el identificador binomial como ultimo nivel a lo largo de la asignacion 
if(identical(bim.obj[,1], taxa.obj[,1])) {
  bim.tax <- data.frame(tax = taxa.obj[,2], binom = bim.obj[,2], stringsAsFactors = FALSE)
  bim.tax <- unite(bim.tax, sep = ";", remove = TRUE, col = Taxonomy ) }
                                              

tax <- strsplit(bim.tax$Taxonomy, ";")
tax <- sapply(tax, "[", c(1:max(lengths(tax)))) # GenBank levels
tax <- as.data.frame(t(tax), stringsAsFactors = FALSE)

# fill NA possition with pervious row assignation;

tax <- data.frame(t(apply(tax, 1, na.locf)))

# also use last_possition_unclassify to fill NA rows <---


#tax <- strsplit(taxa.obj[,2], ";")
#tax <- sapply(tax, "[", c(1:max(lengths(tax)))) # GenBank levels
#tax[is.na(tax)] <- "Unclassified" # fill na possition with tag unknown
#tax <- as.data.frame(t(tax))

dim(tax[complete.cases(tax),]) # 28 with NAs 




Id <- make.unique(as.vector(taxa.obj[,1]), sep = "_")

save <- cbind(Id, unite(tax, sep = ";", remove = TRUE, col = 'Taxonomy'))

save$Taxonomy <- sapply(save$Taxonomy,
                        function(x){gsub(pattern = "$",
                        replacement = ";", x)})


cat('. Formating fasta ...\n')

n_seqs <- length(Id)
seq_headers <- vector(n_seqs, mode="character")

for (i in 1:n_seqs) {
  seq_headers[i] <- paste(">", Id[i], sep = "")
}

fasta <- c(rbind(seq_headers, seqs))

cat('. Writing outputs ... \n')

# 1
write.table(save, file = paste0("Coarbitrator_COI_nuc_curated", ".tax"), append = FALSE, quote = FALSE, sep = " ",
            na = "NA", row.names = FALSE,
            col.names = FALSE)
# 2
write(fasta,
            file=paste0("Coarbitrator_COI_nuc_curated",
                          ".fasta"))

cat('. DONE!\n')

quit(save ='no')
```

 Y verificamos que no tengamos mas duplicados

```bash
awk '{print $1}' Coarbitrator_COI_nuc_curated.tax | sort | uniq -c | sort -n -k2,2 | tail
   # 1 Z93007
   # 1 Z93008
   # 1 Z93009
   # ....
  
grep  -A1 'NC_008833' Coarbitrator_COI_nuc_curated.fasta
>complement(NC_008833
NN
```

error, 

```bash
tail mothur.1557519203.logfile
Z93010 has an error in the taxonomy.  This may be due to a ;;
Z93011 has an error in the taxonomy.  This may be due to a ;;
Z93012 has an error in the taxonomy.  This may be due to a ;;
Z93013 has an error in the taxonomy.  This may be due to a ;;
complement(NC_008833 has an error in the taxonomy.  This may be due to a
```



**use last_possition_unclassify to fill NA rows**

## lineage content

```R
phyla <- data.frame(table(tax$V2))
phyla <- phyla[order(phyla$Freq, decreasing = TRUE), ]
ncol(phyla) # 34


```

Also verify assignation

```R



path_run014 <- '/Users/cigom/metagenomics/COI/co-arbitration-assig/run014_t2_ASVs_99_rdp_outdir'
setwd(path_run014)

fasta.file <- "run014_t2_ASVs.fasta"
tax.file <- "run014_t2_ASVs.Coarbitrator_COI_nuc_curated.wang.taxonomy"


fasta.obj <- readDNAStringSet(fasta.file)
seqs <- as.data.frame(fasta.obj)$x

#   reading FASTA file Coarbitrator_COI_nuc.fasta: ignored 56 invalid one-letter sequence codes

taxa.obj <- read.csv(tax.file, header=FALSE, sep='\t', na.strings=c("","NA"), stringsAsFactors = FALSE)


cat("\n 1. Processing taxonomy... \n")

tax <- strsplit(taxa.obj[,2], ";")
max.rank <- max(lengths(tax))

tax <- sapply(tax, "[", c(1:max(lengths(tax)))) # GenBank levels
tax <- as.data.frame(t(tax))

tax1 <- as.data.frame(apply(tax, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)

rownames(tax1) <- taxa.obj[,1]
                           
cat("\n 2. Bootstrap stats ... \n")

boots <- as.data.frame(apply(tax, 2, function(x) gsub("[A-z||()]", "",  x, perl=TRUE)), stringsAsFactors = F)
boots <- apply(boots, 2, as.numeric)
boots <- data.frame(boots)
boots[is.na(boots)] <- 0
  
rownames(boots) <- taxa.obj[,1]
                             

rank.names <- vector(max.rank, mode="character")

for (i in 1:max.rank) {
  rank.names[i] <- paste("Rank", i, sep="_")
}
                             
colnames(boots) <- rank.names
colnames(tax1) <- rank.names

                             
dim(tax1)
# [1] 21278    20                         
                  
ntaxa <-  function(r) {
  			r <- data.frame(table(r))  
 			  r <- r[order(r$Freq, decreasing = TRUE), ]
  		return(r) }
r1 <- ntaxa(tax1[,1])  

r2 <- ntaxa(tax1[,2]) # 13 - Eukaryota_unclassified
# Eukaryota_unclassified 12, 428 
                             
nrow(tax_clear <- tax1[tax1[,2] != 'Eukaryota_unclassified', ])
# [1] 8850  (41.6 %)

r3 <- ntaxa(tax_clear[,3]) # 31
r4 <- ntaxa(tax_clear[,4]) # 106
r5 <- ntaxa(tax_clear[,5]) #216
r19 <- ntaxa(tax_clear[,19])    

# /Users/cigom/metagenomics/COI/species_resolution_per_db
```

Comparar con asignaciones previas en bold y midori



```bash
/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dada2_asv/run14_20190228_COI/dada2_asv/run014_20190305_COI_t2/run014_t2_ASVs_99_rdp_bold_outdir


```

