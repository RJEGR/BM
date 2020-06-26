Esto es parte del objetivo 

- Evaluación de la técnica en comunidad(es) artificial(es)Definir composición (y función) de mock(s). LE, GD, AS, JR.
- Preparación y secuenciación de mock(s). AS, JR
- Análisis de datos metaB-mock. MM, RG
- wd=`/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish`

# mock metagomfish

La idea es analizar las mock con 2 estrategias.
**1. Datos completos.** Se usan todas las secuencias de la biblioteca. 
**2. Datos filtrados.** Se usan las secuencias filtradas como en [Pruebas parámetros mock](https://app.asana.com/0/827137820949473/1109597644755855). Usando las secuencias de referencia de los organismos contenidos en la mock (`/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/data/run18_20200401_COI`). **Lo que te corresponde hacer [@Ricardo Gore](https://app.asana.com/0/455610172835450/list) es filtrar las bibliotecas.**

![img](https://asana-user-private-us-east-1.s3.amazonaws.com/assets/455610172835444/1117758610266681/bbce67d7a317555e1daa0349f1041bb0?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPD%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQD6zzANxsJ345sw0KYQKhfb9YbJ%2B2%2FKI%2BGrgyRMilDLhQIhAJN8%2Fiko24KTttuf7l7qbQPsBUhwiFlqEZA6y9RYJOyKKrQDCEkQABoMNDAzNDgzNDQ2ODQwIgxsy1TVLYqdFK9QHlcqkQMXuPnf4DkCu3PXBd1XF87%2BZVyq%2F42vYFN2%2Fyk3%2F25zuVM4XvMCtJW8ArbSqPKoSmiAg06gI67Ud1oApnyIXEP%2FgyHuQjEMqs31ONYIPmbMnBJIeSKwoDixleYIZFbf0JDcvSmUKn5pmxzmiZ3HCobXIq%2Bl3CeAw0ZweyZhPbkwLc7dGZfHUMaV1y27pGH9ZqLDhGXA21bvnezs9UAyCeKjuYR5SEYV2N76i0eyW6AVF19DMVqknf%2FAuRD4Qrp4W6sSWstREeAwzf2%2FU%2B%2B3zM8qACP7Cx7W5dnhaVoPkx%2FEGVa13FfnFbMmSxZKMXwwDJqlOoraDfSNPqg93KQo9RozZAwqUUTDzCTCrn0%2BcE3%2BVujKvPGLXbBEUWPuCo6VbNS%2BANviui2P86UevfDPRGwgY1bKV%2F5K58iM753g%2BIM8RNaTAShVoyDxQinYm2Omogyvrqd%2BQoD2mOviJsLEuWfHYRjdCzDjWYA4dRe05r0qehvn7KjBlUnCG4Ah41yWvdy4i5JaptYaqOtwaOo1eOvPwTC745%2F2BTrqAd9%2BDm7uLTIHcUJMtXQ%2F2PzoHndHlfv50HXgK70XjZEkvKZDB%2FVPuTwmCJFQ0VBGLTJT4iV8WZMG23BuF6QVR%2F0k2AL%2F6Flf1y4BbDrXd6wMO7jbji83kpbN9iBcfyaytdQol%2BzWwq%2B2I3HfiwVrdsalokDVa5HoXzZOgLFF5ELn56pN2n1zclqE9R15zbYO7MFE0ncKY0c91Pzl6uo5USlQ81k%2B8q3yi8VjPqTiz04mBiF7YhlhXtW5ZiN8yOfNbeGPZmOP%2BvZ2jYr3vkLy52Xb9y5fSyElz6GuU0iHRDrpeBp97nJkgy2vgQ%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20200522T171330Z&X-Amz-SignedHeaders=host&X-Amz-Expires=120&X-Amz-Credential=ASIAV34L4ZY4OSW532PZ%2F20200522%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=2213fdd9778c1d792aaf3cddf89bb71baa0bb5ff0c97356581f4db0fd25bbdb9#_=_)

1) Modify `*.conf` file

```bash
# This is an example configuration file for FastQ Screen
#

############################
## Bowtie, Bowtie 2 or BWA #
############################
BOWTIE	/LUSTRE/apps/bioinformatica/bowtie1/bowtie
BOWTIE2 /LUSTRE/apps/bioinformatica/bowtie2/bowtie2

############
## Threads #
############

THREADS		24

##############
## DATABASES #
##############

#DATABSE	Relabun ~/metagenomics/metafishgom_mocks/mock_ictio_RA
#DATABS	Congener ~/metagenomics/metafishgom_mocks/mock_ictio_congener
#DATABSE	Diversity ~/metagenomics/metafishgom_mocks/mock_ictio_div

DATABASE	Relabun /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_RA/mock_ictio_RA
DATABASE	Congener /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_congener/mock_ictio_congener
DATABASE	Diversity /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_div/mock_ictio_div
DATABASE Human	/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Genomes/human/hg38ome/hg38
DATABASE Univec	/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Genomes/UniVec/UniVec
```

2) Check the quality of raw sequence

```bash
mkdir fastqc

fastqc *.gz -t 24 --nogroup -o ./fastqc &

mkdir multiqc
export PATH=/LUSTRE/apps/Anaconda/conda2/bin:$PATH
source activate multiqc_py2.7

multiqc ./fastqc/*zip -o ./multiqc
```

3) Then, set the flags:

```bash
# Para esta prueba usamos el algoritmo bowtie2 atraves de fastqscreen 

# Consider mapping to three genomes (A, B and C), the string '003' produces a file in which reads do not map to genomes A or B, but map (once or more) to genome C.  The string '--1' would generate a file in which reads uniquely map to genome C. Whether reads  map to genome A or B would be ignored.

# 1. Fix headers name due to fastq_screen 

for i in $(ls *R1*fastq.gz); do zcat $i | awk '/^@M/{print $1; next}{print}' > ${i%-*}_R1.fastq ; done

# Also for Reverse reads
for i in $(ls *R2*fastq.gz); do zcat $i | awk '/^@M/{print $1; next}{print}' > ${i%-*}_R2.fastq ; done

```

Then,

```bash

# Sanity check, reads must be equals!

for i in $(ls *.fastq.gz); do zcat $i | grep -c '^@M'; done
for i in $(ls *.fastq); do cat $i | grep -c '^@M' $i; done

```

Next step

```bash
#!/bin/bash
#
#SBATCH -J fscreen
#SBATCH -p cicese
#SBATCH -t 6-00:00:00
#SBATCH -N 1
#SBATCH -n 24

conf=$1
outdir=$2
fastq=$3

mkdir -p $outdir

fastq_screen --aligner bowtie2 --conf $conf --subset 0 --filter 3-- --tag --force --outdir $outdir $fastq

# 3. Fix headers name due to fastq_screen

filtered=$outdir/${fastq%.fastq}*.tagged_filter.fastq
hits=`basename $filtered`
mock_hits=${hits%.tagged_filter.fastq}.mock_hits.fastq

awk '/^@M/{gsub(/#FQST%*/, " "); print ""$1; next}{print}' < $filtered > $outdir/$mock_hits

# 4.  Parse R2 reads based on R1 filtered reads 

awk '/^@M/{print $1}' < $outdir/$mock_hits > ${mock_hits%.fastq}.ids

# 5. get sequences based on a list

file_R2=${fastq%_R1.fastq}_R2.fastq
list=${mock_hits%.fastq}.ids

grep -A 3 -Ff $list $file_R2 | sed '/^--$/d' > $outdir/${file_R2%.fastq}.mock_hits.fastq

awk '/^@M/{print $1}' < $outdir/${file_R2%.fastq}.mock_hits.fastq > ${file_R2%.fastq}.mock_hits.ids

diff *ids

rm *ids

exit

# for i in $(ls *Con*R1*); do sbatch fscreen.slurm congener.conf mock_ictio_congener $i; done
# for i in $(ls *RelA*R1*); do sbatch fscreen.slurm RA.conf mock_ictio_RA $i; done
# for i in $(ls *Div*R1*); do sbatch fscreen.slurm div.conf mock_ictio_div $i; done



```

Further

```bash
# Additionals
# get back the hit_no_genome reads
fastq_screen --conf ../fastqscreen.conf --nohits 012-X04-P68-COI-ICT_S46_L001_R1_001.fastq.headers.tagged.fastq
no_hits=012-X04-P68-COI-ICT_S46_L001_R1_001.fastq.headers.tagged.tagged_filter.fastq

awk '{s++}END{print s/4}' $no_hits

awk '/^@M/{gsub(/[#FQST]/, " "); print ""$1; next}{print}' < $no_hits > 012-X04-P68-COI-ICT_S46_L001_R1_001.no_hits.fastq
no_hits=012-X04-P68-COI-ICT_S46_L001_R1_001.no_hits.fastq

grep "^@M" -A 1 $no_hits | sed '/^--$/d' | sed 's/^@/>/g' > 012-X04-P68-COI-ICT_S46_L001_R1_001.no_hits.fasta

# Extra. get Reverse sequences  from no_genome reads

no_hits_R1=012-X04-P68-COI-ICT_S46_L001_R1_001.no_hits.fastq

awk '/^@M/{print $1}' < $no_hits_R1 > 012-X04-P68-COI-ICT_S46_L001_R1_001.no_hits.ids

list=012-X04-P68-COI-ICT_S46_L001_R1_001.no_hits.ids


grep -A 3 -Ff $list $file_R2 | sed '/^--$/d' > 012-X04-P68-COI-ICT_S46_L001_R2_001.no_hits.fastq

awk '/^@M/{print $1}' < 012-X04-P68-COI-ICT_S46_L001_R2_001.no_hits.fastq > 012-X04-P68-COI-ICT_S46_L001_R2_001.no_hits.ids

diff *no_hits.ids
```



Archivos:

(P68) Drive-slides **Mock usage**

(P68) run12_P68_resultados

~/cigom/paper_metafishgom/mock/**mock_coi_peces.log** (solo peces)

~/cigom/log_and_notes/mocks_coi.log (todo menos peces)

[por hacer] metafishgom_disponible_para_mock.tsv

Archivo con los linajes worms de los organismos disponibles para hacer mocks

Esta parte se complementa con [Primer mismatches](https://app.asana.com/0/0/1152527708347739):

```bash
Al archivo de targetDB. Aniadir una columna por primer. Donde se indique cuantos mismatches tiene la secuencia con el primer. Usar NA en caso que la region no este cubierta.

Leray_2013 dice que encontraon maximo 6 mismatches analizando organismos de 6 filos. Ver Figura 3.

[Revisar si existe alguna relacion entre grupos biologicos y numero de mismatches]

Alternativa 1. Revisar primer-mismatch solamente en las secuencias de los organismos que se tienen disponibles para hacer la mock de interes. Generar un archivo con los linajes (ver mock metafishgom)
Alternativa 2. Usar mitogenomas (que tienen el COI completo) para hacer el calculo de primer-mismatch por grupo

Objetivos:
1. Cuantificar el numero de mismatches por grupo taxonomico
2. Evaluar primer efficiency (bias)

Meta: Se hara una mock que contenga especies con diferente grado de mismatch (usando misma cantidad de DNA) para evaluar si hay sesgo en la amplificacion.

Opcion propuesta es DECIPHER::designprimers y relacionados
```

# mock_congener

Arreglo de organismos en Drive cigom_ICT_12S > mock_congener

https://docs.google.com/spreadsheets/d/1gbed0kfpDwu1pNL0SsX_XimXavqD4WVlyf9T4xq

[6kUY/edit?usp=sharing](https://docs.google.com/spreadsheets/d/1gbed0kfpDwu1pNL0SsX_XimXavqD4WVlyf9T4xq6kUY/edit?usp=sharing)

Marcador 12S secuenciado en [✓ Analisis de run17 con repo mothur_18S](https://app.asana.com/0/0/1166342506576482)



Refseqs:

/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/procesamiento_ecologico/metafishgom_mocks/mock_ictio_congener_refseqs.fasta

# mock relabun

Refseqs:

/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/procesamiento_ecologico/metafishgom_mocks/mock_ictio_RA_refseqs.fasta



# mock diversidad




Refseqs:

/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/procesamiento_ecologico/metafishgom_mocks/mock_ictio_RA_refseqs.fasta





In R

```R
library(tidyverse)
library(ggplot2)

pct <- read.csv('fastqscreen_pct.csv', sep=',', stringsAsFactor =F)
reads <- read.csv('fastqscreen_reads.csv', sep=',',stringsAsFactor =F)

pct %>% pivot_longer(cols = contains("hit")) %>%
filter(!Genome %in% c('Human', 'Univec')) %>%
ggplot(aes(x = Genome, y = Unmapped, fill = name)) +
  geom_bar(stat = "identity", position = "stack", color = 'black') +
	facet_grid(~Rep, scales = 'free', space='free') +
	labs(x = 'Mock', y = 'Reads (Pct)') +
	scale_fill_grey() +
	theme(axis.text.x = element_text(angle = -45, hjust = 0), 
              axis.text.y = element_text(size  = 7))
```



## Check Coverage

```bash
# To assess the read composition of our assembly, we want to capture and count all reads that map to our assembled transcripts, including the properly paired and those that are not we run the process below.

#DATABASE	Relabun 
/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_RA/mock_ictio_RA
#DATABASE	Congener 
/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_congener/mock_ictio_congener

#DATABASE	Diversity 
/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_div/mock_ictio_div

#!/bin/bash
#
#SBATCH -J bowtie
#SBATCH -p cicese
#SBATCH -t 6-00:00:00
#SBATCH -N 1
#SBATCH -n 24

index=$1 # refseq indexed
forward=$2 # 018-Mock-Con1-COI_R1.mock_hits.fastq
reverse=$3 # 018-Mock-Con1-COI_R2.mock_hits.fastq
outTag=${forward%_R1.mock_hits.fastq}

# 2
# suppress SAM records for unaligned reads --no-unal
# srun bowtie2 -p 12 -q --no-unal -k 20 -N -x $index -1 $forward -2 $reverse | samtools view -@12 -Sb -o ${outTag}.mock_hits.bam

bowtie2 -p 12 -q --no-unal -k 20 -N 1 --no-discordant -x $index -1 $forward -2 $reverse | samtools view -@12 -Sb | samtools sort -o ${outTag}.mock_hits.sorted.bam

# remove dup
# samtools rmdup -S no_discordant.sorted.bam no_discordant.sorted.rmdup.bam


exit
# 3

#bamtools stats -in ${outTag}.mock_hits.bam

# Use in the loops
for i in $(ls *Con*R1*); do sbatch reads_represented.slurm /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_congener/mock_ictio_congener $i ${i%_R1.mock_hits.fastq}_R2.mock_hits.fastq; done
#
for i in $(ls *RelA*R1*); do sbatch reads_represented.slurm /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_RA/mock_ictio_RA $i ${i%_R1.mock_hits.fastq}_R2.mock_hits.fastq; done
#
for i in $(ls *Div*R1*); do sbatch reads_represented.slurm /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/mock_metagomfish/db2screen/mock_ictio_div/mock_ictio_div $i ${i%_R1.mock_hits.fastq}_R2.mock_hits.fastq; done


# bam2sam
# samtools view -h -o 018-Mock-Con1-COI.mock_hits.sam 018-Mock-Con1-COI.mock_hits.bam
# for i in *bam; do samtools view -h -o ${i%.bam}.sam $i ; done
# bam2sam (only uniques)
# for i in *bam; do samtools rmdup -S $i -| samtools view -h -o ${i%.bam}.rmdup.sam; done

#

for i in *rmdup.sam; do awk '{print $3}' $i | grep -v -e '^LN:' -e '^SO:' -e '^PN:' | uniq -c ; done

# indexes
# samtools index 018-Mock-Con1-COI.mock_hits.sorted.bam
for i in *sorted.bam; do samtools index $i; done
```

Then, in R

```R
BiocManager::install("ggbio")
library(ggbio)
library(GenomicRanges)
library(Biostrings)

fasta <- dir(pattern='fasta', path='~/metagenomics/metafishgom_mocks', full.names=TRUE)

ref <- readDNAStringSet(fasta[1])
widths <- width(ref)
GC <- rowSums(letterFrequency(ref, letters="GC", OR=0)) / widths

id <- sapply(strsplit(names(ref), " "), `[`, 1)

wh <- GRanges(seqnames = Rle(id),
    ranges = IRanges(1, end = widths, names = id),
    Rle(strand(c("+"))),     
    GC = GC)

# compute the GC content in a sliding window (as a fraction) for a sequence no. 364
# window = 100
#gc <- rowSums(letterFrequencyInSlidingView(ref[[1]], window, c("G", "C")))/window
#plot(gc, type = 'l')

fl.bam <- dir(pattern="bam", path='~/metagenomics/metafishgom_mocks/coverage/bam', full.names=TRUE)

options(stringAsFactor = FALSE)

# autoplot(fl.bam[1], which = wh)

library(Rsamtools)

bam <- scanBam(fl.bam[1])
names(bam[[1]]) # 

#function for collapsing the list of lists into a single list

.unlist <- function (x){
   ## do.call(c, ...) coerces factor to integer, which is undesired
   x1 <- x[[1L]]
   if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
   } else {
      do.call(c, x)
   }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

               
#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)
# refseq coverage
sum(id %in% names(table(bam_df$rname)))
#how many entries on the negative strand of chr22?
rname <- 'ConsensoIctio65'
table(bam_df$rname == rname & bam_df$flag == 16)
# FALSE    TRUE 
#3875997   24413

#function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

#test neg function with subset of chr22
test <- filter(bam_df, rname == rname)
dim(test)
#[1] 56426    13
table(apply(as.data.frame(test$flag), 1, check_neg))
#number same as above
#FALSE  TRUE 
#32013 24413

#function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}

#check pos function
table(apply(as.data.frame(test$flag), 1, check_pos))
#looks OK
#FALSE  TRUE 
#24413 32013

#store the mapped positions on the plus and minus strands
rname_neg <- bam_df[bam_df$rname == rname &
                    apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos'
                   ]
length(rname_neg)
#[1] 24413
rname_pos <- bam_df[bam_df$rname == rname & apply(as.data.frame(bam_df$flag), 1, check_pos), 'pos' ]
length(rname_pos)
#[1] 32013
 
#calculate the densities
rname_neg_density <- density(rname_neg)
rname_pos_density <- density(rname_pos)
 
#display the negative strand with negative values
rname_neg_density$y <- rname_neg_density$y * -1
 
plot(rname_pos_density,
     ylim = range(c(rname_neg_density$y, rname_pos_density$y)),
     main = "Coverage plot of mapped CAGE reads",
     xlab = rname,
     col = 'blue',
     type='h')
lines(rname_neg_density, lwd=2.5, col = 'red', type='h')
               
legend(2,1,c("+ strand","- strand"), lwd=c(5,2), col=c("blue","red"), y.intersp=1.5)

```



## Remove duplicates

Consideramos que cada lectura con puntaje de calidad phred >= 20 es una secuencia con variante única. Usamos esta etapa para remover duplicados en pruebas que permitan usar una sola direccion de la lecturas paired-end. Por ejemplo:

```bash
# https://bioinf.shenwei.me/seqkit/usage/
seqkit rmdup -s 018-Mock-Con1-COI_R1.mock_hits.fastq > 018-Mock-Con1-COI_R1.mock_hits.rmdup.fastq # 14607 duplicated records removed
seqkit rmdup -s 018-Mock-Con1-COI_R2.mock_hits.fastq > 018-Mock-Con1-COI_R2.mock_hits.rmdup.fastq # 13483 duplicated records removed
```

