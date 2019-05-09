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

