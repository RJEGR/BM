## Summary

En este caso (Febrero 2020) vamos a construir el modelo usando como referencia el fragmento de Pez Cebra (`COI_danio_rerio.fasta`) y como secuencias a alinear una lista de  3825 mitogenomas descargados de ncbi en vez de centroides del set de NCBI.

sobre los mitogenomas realizasamos la siguiente limpieza:

- Removemos abreviaciones raras
  - 144 / 3825 secuencias con abreviaciones raras
  - retenemos la lista con: `seqkit`
- Dereplicamos con search
  - 2102 / 3684 unique sequences
- Revision del alineamiento del modelo obtenido con `mafft`: (Remover inserciones intergenicas (`NNN` y _gaps_)
  -  1768 / 2103

Los archivos de entrada se encuentran en:

```bash
# workdir: /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/ncbi/COI_de_mitogenomas_clean
# 1 /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/coi_refaln/

f12a0653af018d7d839d7d2b68e17acd  COI_danio_rerio.fasta
d2548476f3313fdfc4e9c4323c0046a2  lscoi_drerio_ictioconsenso_primers_aln.fasta
4ce1c65c46bbb0e21ddddc46ab9a19ac  COI_danio_rerio_coords.txt

# 2 /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/ncbi/COI_de_mitogenomas

68fa9f7cf4cd0bf22efb398a8a6df35d  ncbi_fish_mitogenomes_Feb2020.fasta
```

0. Removemos abreviaciones raras

```bash
fasta=ncbi_fish_mitogenomes_Feb2020.fasta


grep '^>' $fasta | sed 's/>//g' > ${fasta%.fasta}.headers.tmp

egrep 'sp\.|cf\.|aff\.|hybrid|UNVERIFIED|environmental|[a-z][[:blank:]]x[[:blank:]][A-Z]' ${fasta%.fasta}.headers.tmp | cut -d' ' -f1 > ${fasta%.fasta}.redundant.tmp

egrep -v "sp\.|cf\.|aff\.|hybrid|UNVERIFIED|environmental|[a-z][[:blank:]]x[[:blank:]][A-Z]" ${fasta%.fasta}.headers.tmp | cut -d' ' -f1 > ${fasta%.fasta}.good.tmp

wc -l *tmp

# in cluster select good files

cat $fasta | seqkit grep -f ${fasta%.fasta}.good.tmp > ${fasta%.fasta}.good.fasta

# clean
rm *tmp
```

1. Dereplicamos

```bash
seqs=${fasta%.fasta}.good.fasta
sorted=${seqs%.fasta}_sorted.fasta
derep=${seqs%.fasta}_derep.fasta

vsearch  -sortbylength $seqs -output $sorted
vsearch -derep_fulllength $sorted -output $derep -sizeout
```

2. Alineamos todas las secuencias: `sbatch mafft.sh ncbi_fish_mitogenomes_Feb2020.good_derep.fasta COI_danio_rerio.fasta`

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

**OPTIONAL:** Using Decipher to construc model based on the Open Reading Frame

```bash
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

