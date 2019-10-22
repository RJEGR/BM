# NCBI Entrez Direct UNIX E-utilities 

It is available to download from the NCBI website [here](ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect) or [here](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

Dowload full ncbi COI barcodes data-base using ncbi nucleotide search as Teresita et.al (2018):

> CO1[GENE] OR COI[GENE] OR COX1[GENE] OR COXI[GENE] AND "Eukaryota"[Organism] 

**La descarga se realizo con fecha Octuber 08, 2019:**

- Downolad method: **e-utilities**

- Al usar el tag "BARCODE" en la busqueda de secuencias despues del 2015: **100,273**
- Al filtrar secuencias con abreviaciones redundante en nivel especie: **89,889**
- Al filtrar secuencias menores o igual a 200 pares de bases: **89,847**
  - Se cuenta con informacion taxonomica para estas secuencias: **VERDADERO**
- Al filtrar mitogenomas (secuencias mayores a 2000 pb) : **88,394**
- al de-replicar secuencias: **54,171**

## Add API to faster access

```bash
NCBI_API_KEY=74e0e4bf2d8eaa3bb742c46316dbafe12909
echo "NCBI_API_KEY=74e0e4bf2d8eaa3bb742c46316dbafe12909" >> $HOME/.bash_profile
```

## Descarga de secuencias

```bash
esearch -db nucleotide -query "CO1[GENE] OR COI[GENE] OR COX1[GENE] OR COXI[GENE] OR BARCODE AND "txid7898"[Organism:exp] AND COI AND "2016/01/01:2019/09"[MDAT]" | efetch -format fasta > CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta
```

## Removiendo Abreviaciones redundantes:

Taxa were filtered according to the concent of the species field so that only fully identified taxa with a complete latin binomial (genus and species) were retained. Entries that contained the abbreviations. `sp.`, `nr.`, `aff.`, or `cf.` were discarded (Teresita porter et al 2018). Also include abbreviations `UNVERIFIED` and `environmental` (Miguel M 2019).

1) Wrangling:

```bash
# 1)
grep '^>' BARCODE_txid7898_CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta | sed 's/>//g' | awk '{print $2,$3, $4}' | sort | uniq -c | sort -n > frequency.txt

# number of redundant abbreviations: 10361

egrep 'sp\.|cf\.|aff\.|UNVERIFIED|environmental' frequency.txt | awk '{ sum += $1; } END { print sum; }' "$@"

# Linearize fasta to grep inVerte and format back to fasta
awk -f linearizefasta.awk < BARCODE_txid7898_CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta | egrep -v 'sp\.|cf\.|aff\.|UNVERIFIED|environmental' | tr "\t" "\n" > ncbi_complete.fasta

# Sanity check
grep -c "^>" ncbi_complete.fasta # 89889 

```

2.1) seqmagick —min-length 200 max length

```bash
seqmagick mogrify --min-length 200 --max-length 2000 ncbi_complete.fasta
grep -c "^>" ncbi_complete.fasta # 89847

# Sanity check
egrep -c 'sp\.|cf\.|aff\.|UNVERIFIED|environmental' ncbi_complete.fasta
```

Save headers

```bash
grep '^>' ncbi_complete.fasta | sed 's/>//g' > ncbi_complete.headers
```

0) Additional 

```bash

freq=frequency.txt
awk '{ sum += $1; } END { print sum; }' "$@" $freq # 100273

# 7026
grep 'sp\.' $freq | awk '{ sum += $1; }END{ print sum; }' "$@" 

# 1044
grep 'cf\.' $freq | awk '{ sum += $1; }END{ print sum; }' "$@" 

# 191
grep 'aff\.' $freq | awk '{ sum += $1; }END{ print sum; }' "$@" 

# 1363 <- incluye sp
grep 'UNVERIFIED' $freq | awk '{ sum += $1; }END{ print sum; }' "$@" 

# 1237
grep 'environmental'  $freq | awk '{ sum += $1; }END{ print sum; }' "$@" 

grep 'nr\.' $freq # 0

# 7026 + 1044 + 191 + 1363 + 1237 + 0 = 10861 <- here we have duplicate 

```

Save incomplete

```bash
awk -f linearizefasta.awk < BARCODE_txid7898_CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta | egrep 'sp\.|cf\.|aff\.|UNVERIFIED|environmental' | tr "\t" "\n" > ncbi_incomplete.fasta

grep -c "^>" ncbi_incomplete.fasta # 10384

grep '^>' ncbi_complete.fasta | sed 's/>//g' | awk '{print $2,$3}' | sort | uniq -c | sort -n > frequency_sp.txt

```

`linearizefasta.awk`

```bash
/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;}
     {printf("%s",$0);}
END  {printf("\n");}
```



## De-replicate sequences

```bash
vsearch -derep_fulllength ncbi_complete.fasta -output ncbi_complete_drep.fasta
grep -c "^>" ncbi_complete_drep.fasta
```

## Primer design

## Aligment

1)

```R
# load the DECIPHER library in R
library(DECIPHER)

getwd()

# specify the path to the FASTA file (in quotes)
fas <- "ncbi_complete.fasta"

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet(fas)

# look at some of the sequences (optional)
seqs

# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
# seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned,
   file="ncbi_complete.align")

quit(save = 'no')
```

2)  muscle For low number of sequence (`Segmentation fault error if thousan of sequences`)

```bash
srun muscle -in ncbi_complete.fasta -out ncbi_complete.muscle.align &
```

3) mafft

```bash
#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

input=$1
output=${input%.fasta}.mafft.fasta

mafft --thread $SLURM_NPROCS $input > $output

exit
```

## Get the complete taxonomy

Esta es una version para recuperar el linaje de alguna secuencia del ncbi basado en el identificador gi y usando la agrupacion de generos.

```bash
awk '{print $1}' ncbi_complete.headers | sed 's/>//'  > ids

# or separate by genus taxons
cat ncbi_complete.headers | awk '{gsub(/[[:punct:]]/, "" , $2); print $2}' | sort | uniq -c | sort -k2,2 > ncbi_complete_genus

wc -l ncbi_complete_genus # ~ 3192 ncbi_complete_genus
```

Then grep genus by genus

```bash
cat ncbi_complete_genus | awk '{print $2}' > genusList
# Debido al procesamiento de usar la segunda columna para agrupar generos, perdemos detalles del identificador acc MH925110.1, por lo que buscamos manualmente su linaje y adjutamos al final dentro del archivo >> accTaxId
# Make the ids files per genus . time demand

while IFS= read -r pattern; do grep $pattern ncbi_complete.headers | awk '{print $1}' > ${pattern}.ids & done < genusList; wait; echo "Grouping  genus done"

# Sanity check 
ls *ids | wc -l # ~ 3192 genus files
wc -l *ids | sort -k2,2 | sed 's/.ids//g' > genusFiles

diff genusFiles ncbi_complete_genus  # must be zero
# tenemos un problema con el reconoimiento de patron 
# time demand
for i in *ids
do 
	genus=${i%.ids}
	epost -db nuccore -format acc -input $i | \
	esummary | \
	xtract -pattern DocumentSummary -element AccessionVersion TaxId > ${genus}.taxId
done

# Remove empty files
find . -size 0 -delete

# Sanity check

ls *taxId | wc -l # must be ~ 3193
wc -l *taxId | sort -n | head

# checkpoint
ls *taxId | awk '{gsub(/[.]/, " " , $0); print $1}' | sort > DowloadedList

diff DowloadedList genusList | awk '">" {gsub(/>/,"", $0); print}' | sort > checkpointList

wc -l checkpointList # must be equal zero

for genus in $(cat checkpointList)
do 
	epost -db nuccore -format acc -input ${genus}.ids | \
	esummary | \
	xtract -pattern DocumentSummary -element AccessionVersion TaxId > ${genus}.taxId
done

# redoing sanity check until checkpointList must be zero
# Then,
# get taxonomy

cat *taxId | awk '{print $2}' | sort | uniq > completeTaxIds # 10866 number of genus

for i in $(cat completeTaxIds); do efetch -db taxonomy -id $i -format xml | \
xtract -pattern Taxon -tab "," -first TaxId ScientificName \
-group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
-block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
-block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
-block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
-block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
-block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
-block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
-group Taxon -tab ";" -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS"; done > completeTaxIds.taxonomy

# Sanity check
wc -l completeTaxIds.taxonomy # 10866 number of linage recovery
cut -d',' -f1 completeTaxIds.taxonomy | sort > DowmloadedLinage
diff DowmloadedLinage completeTaxIds | sort # must be zero

# Check Actinopteri class
grep 'Actinopteri' -c completeTaxIds.taxonomy # 10862

cut -d',' -f5 completeTaxIds.taxonomy | sort | uniq -c | sort -n
#   1 Chordata
#   4 Cladistia
#10861 Actinopteri

# clean temporally files
rm *.ids

cat *.taxId | sort | uniq > accTaxId


```

Encontraremos el error de desfase de los linajes, por ejemplo tenemos 579 ordenes con el símbolo `-` en el nivel orden. 

Parse taxonomy to fasta and separe fasta by group

```bash
awk -f linearizefasta.awk < ncbi_complete.fasta | sort -k1,1 | sed 's/^>//g' > ncbi_complete.fasta.table

# accTaxId
# parse both files:

# accTaxId ----
# AB042837.1	7904

# completeTaxIds.taxonomy ----
# 7904 Acipenser transmontanus,Metazoa,Chordata,Actinopteri,Acipenseriformes,Acipenseridae,Acipenser

# Final output ----

# AB042837.1	7904 Acipenser transmontanus,Metazoa,Chordata,Actinopteri,Acipenseriformes,Acipenseridae,Acipenser

```

## Split taxonomy to sequence data in R

```R
library(dplyr)
x <- read.table('accTaxId', header = FALSE)
y <- read.table('completeTaxIds.taxonomy', header = FALSE, sep=',')

names(x) <- c('acc','taxid')
names(y) <- c('taxid', 'SPECIE','KING', 'PHYL', 'CLSS', 'ORDR', 'FMLY', 'GNUS')

x %>% inner_join(y) %>% as_tibble() %>% select(acc, taxid, KING, PHYL, CLSS, ORDR, FMLY, GNUS, SPECIE) -> tax

library(Biostrings)

z <- readDNAStringSet('ncbi_complete.fasta')
names <- sapply(strsplit(names(z), " "), `[`, 1)
names(z) <- names

z <- as.data.frame(z[order(names(z))])
z$acc <- rownames(z)

# processing fasta to save  

tax %>% inner_join(z) %>% as_tibble() -> save
save$acc <- sub("^", ">", save$acc)

# unite(1:9, col='id', sep = "|")

# save fasta
wf <- unite(save, 1:9, col='id', sep = "|")
wf <- c(rbind(wf$id, wf$x))
write(wf, file=paste0("ncbi_complete.fasta.bkp"))

# Only sequence

ws <- c(rbind(save$acc, save$x))
write(ws, file=paste0("ncbi_complete.fasta"))

# Only taxonomy
wt <- unite(save, 3:9, col='tax', sep = ";") %>% select(acc, tax) %>% as.data.frame()
wt$acc <- sub("^>", "", wt$acc)
wt$tax <- sub("$", ";", wt$tax)
#names(wt) <- NULL
write.table(wt, file=paste0("ncbi_complete.taxonomy"), sep=" ", 
            row.names = FALSE, 
            col.names = FALSE,
            quote=FALSE)

# Plot a visualization

library(data.table)
taxtb <- data.table(table(tax$GNUS))
names(taxtb) <- c("rank", "n")
# cumulative distribution of reads
par(mfrow=c(1,2))
plot(ecdf(taxtb$n), main = "Cumulative distribution of genus", xlab='Number of sequences')
sample <- filter(taxtb, n <= 100)
plot(ecdf(sample$n), main = "Cumulative distribution above the threshold", xlab='Number of sequences')
#taxtb$n <- taxtb$n / sum(taxtb$n) *100
taxtb <- data.table(table(tax$ORDR))
names(taxtb) <- c("rank", "n")
taxtb <- taxtb[order(-n), ]
taxtb[n <= 100, rank := "Others"]

library(ggpubr)

ggbarplot(taxtb, x = "rank", y = "n",
#          palette = "Paired",            # jco journal color palett. see ?ggpar
          x.text.angle = 90,           # Rotate vertically x axis texts
          ylab = "Number of sequence",
          xlab = "Orders",
          rotate = TRUE,
          ggtheme = theme_minimal()
          #facet.by = "Type"
          ) + theme(axis.text.y = element_text(hjust = 1, size = 7))


# ggalluvial (time demand) 
library(ggalluvial)
alluv <- to_lodes_form(data.frame(tax), key = "Rank", axes = 4:7)
# alluv <- filter(alluv, Rank == 'ORDR')

ggplot(data = alluv,
       aes(x = Rank, stratum = stratum, alluvium = alluvium,
           label = stratum)) +
  geom_stratum() + 
  geom_text(stat = "stratum", size = 3) +
  geom_flow(stat = "alluvium",
            aes.bind = TRUE, lode.guidance = "rightward") +
  theme_minimal() 
#+
#  ggtitle("") +
#  xlab("") + ylab("")

```



### Retrieve Taxon IDs from list of genome accession number

> From https://github.com/NCBI-Hackathons/EDirectCookbook

```bash
# 0) 
awk '{print $1}' ncbi_complete.headers | sed 's/>//' > ids
#

# 1 get the taxid per species
cat ids | \
epost -db nuccore -format acc | \
esummary | \
xtract -pattern DocumentSummary -element AccessionVersion TaxId \
> taxId
```

```bash
for i in $(awk '{print $2}' Astyanax.taxId); do efetch -db taxonomy -id $i -format xml | \
xtract -pattern Taxon -tab "," -first TaxId ScientificName \
-group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
-block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
-block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
-block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
-block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
-block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
-block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
-group Taxon -tab "," -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS"; done
```

## First intention (Above agoust 2019)

- Download Date: July/9/2019

- Number of sequence: 3,388,817 (2,099,002 uniques)
- Size : **5.2G** 
- Filename: CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta
- **Cluster path**: /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/ncbi/CO1_COI_COX1_COXI_GENE_Eukaryota_nr/CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta
- Downolad method: **e-utilities



## Align database

In order to compare the complete folmer and leray primer intraspecific and interspecific and genetic distance by groups of zooplankton we set and aligment of the database.

1) First lets make uniques barcodes in the database:

```bash
#!/bin/bash
#SBATCH --job-name=uniques
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=12

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

# Variables

fasta=$1 # sequence to align

mothur="/LUSTRE/bioinformatica_data/genomica_funcional/bin/mothur/mothur"

$mothur "#summary.seqs(fasta=$fasta);unique.seqs(fasta=current);summary.seqs(fasta=current)"

exit
```

2) Then, lets make a pair-wise alignment of the unique sequences:

```bash
#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"

./clustalo -i CO1_COI_COX1_COXI_GENE_Eukaryota_nr.unique.fasta -o CO1_COI_COX1_COXI_GENE_Eukaryota_nr.unique.align --threads $SLURM_NPROCS

exit

```

3) Cut the leray fragment

> primer 
>
> GGWACWGGWTGAACWGTWTAYCCYCC 
>
> TAIACYTCIGGRTGICCRAARAAYCA

```bash
#!/bin/bash
#SBATCH --job-name=pcr.seqs
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

fasta=$1 # CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta

mothur="/LUSTRE/bioinformatica_data/genomica_funcional/bin/mothur/mothur"

$mothur "#pcr.seqs(fasta=$fasta, oligos=coi.oligos, processors=$SLURM_NPROCS);summary.seqs(fasta=current)"

exit

```

3.1) Finaly, lets cut the leray fragment (F/R primers)

```bash
#!/bin/bash
#SBATCH --job-name=pcr.seqs
#SBATCH -N 3
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

# Variables

fasta=$1 # CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta
ref=$2 # CO1_COI_COX1_COXI_GENE_Eukaryota_nr.unique.align

mothur="/LUSTRE/bioinformatica_data/genomica_funcional/bin/mothur/mothur"

$mothur "#align.seqs(fasta=$fasta, reference=$ref, processors=$SLURM_NPROCS);summary.seqs(fasta=current)"

exit
```



### Notes

*Even when using GenBank, which is the most redundant sequence repository, the coverage of target species is not satisfactory, as the available reference sequences cover only a fraction of marine zooplankton ([sergio stefanni et. al 2017](https://www.nature.com/articles/s41598-018-30157-7)).*