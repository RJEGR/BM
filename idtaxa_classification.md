El algoritmo esta divido en dos fases: 

- *Learn-taxa* 
- *IdTaxa functions*



El primer modulo (learn-taxa) trabaja con secuencias asigandas y corrige falsas asignaciones con una estrategia denominada *tree-descent* (una aproximacion al arbol de desiciones del Machine learning) agrupando k-meros de secuencias en una raiz de asignacione en comun (ver diagrama 1S del paper), a esto le llaman los autores la etapa de entrenamiento, posteriormente, se alimenta otro grupo de secuencias sin clasificacion para la etapa de asignaciones (IDtaxa) en donde utiliza la informacion de la etapa de entrenamiento para asignar.

> El proposito de la etapa learn-taxa es identificar errores en la taxonomia (training set) y ascelerar la clasificación de las secuencias *query* con la funcion iDtaxa.
>
> La funcion *idTaxa* toma como entrada el objeto de salida de la etapa learn-taxa asi como las secuencias *query* para clasificar. Esta funcion tiene como salida la clasificacion de cada secuencia *query* en el formato de asigacion taxonómica con su valor de confiabilidad para cada nivel-rank.

El trabajo de entrenamiento de la base de datos BOLD tomo mucho tiempo y aborto desde un equipo mac (16GB de ram, 3.1 GHz Intel Core i7); es necesario probar en el cluster una vez se definan los parametros optimos (threshold, maxGroupSize,  maxIterations, allowGroupRemoval).

Se puede pre-procesar las bases de datos (ej. midori, bold) con la fase de aprendizaje *learnTaxa* de manera que se identifican errores en la taxonomia y optimiza el tiempo de corrida de la asignacion.



Paso 1:

Dentro del cluster ejecutamos la etapa `learnTaxa`; algunos parametros a configurar pueden ser los enlizatados a continuacion, sin embargo se usaron los default:

1. maxGroupSize, 
2. maxIterations, 
3. allowGroupRemoval

Es necesario ademas, establecer el rank Root, en el archivo de taxonomía.  En el caso de esta prueba, se uso una base de BOLD previamente procesada y fue necesario hacer la siguiente modificación: `sed 's/root;/Root;/g' BOLD_public_species.tax > BOLD_public_species.tax.Root`

Posteriormente corremos el siguiente script en R desde el modulo:

`module load R-3.5.0`

```R
rm(list=ls());

# if RData is in directory continue with LearnTaxa step

# Path ----
library(DECIPHER)

path = getwd()

setwd(path)

##  args insertion ----

# args = commandArgs(trailingOnly=TRUE)

#if (length(args)<3) {
#  stop("!!!\n. 
#       PLEASE, INPUT NECESSARY FILES IN THE SYNTAXIS AS FOLLOW EXAMPLE:\n
#       Rscript --vanilla IDTaxa.R input[string] reference[string] taxonomy[string] #threshold[interger] .\n", call.=FALSE)
#} else { 
#  query <- args[1]
#  seqs_name <- args[2]
#  tax_name <- args[3] }

query <- "OMEGA_A_1e_120_asvs.fasta"
seqs_name <- 'BOLD_public_species.fasta'
tax_name <- 'BOLD_public_species.tax.Root'

# Read inputs ----

q_path <- paste0(path, "/", query)

# Sequnces references path
seqs_path <- paste0(path, "/", seqs_name)

# taxonomy assigments path
tax_path <- paste0(path, "/", tax_name)

#
tag <- strsplit(query, "[.]")[[1]]
tag_ref <- strsplit(seqs_name, "[.]")[[1]]
out_name <- paste0(paste(tag[-length(tag)], collapse="."), "_",
                   paste(tag_ref[-length(tag_ref)], collapse="."), "_",
                                                  "DECIPHER.taxonomy")
#
# Parameters ----
# 
# if (length(args) > 3) { threshold <- as.numeric(args[4]) } else { threshold <- 99 }


## Training parameters

maxGroupSize <- 1 # max sequences per label (>= 1); can be set to Inf (infinity) to allow for an unlimited number of sequences per group.
maxIterations <- 1 # must be >= 1; also can remove mislabeled in the training data by turning itinerations > 1 
allowGroupRemoval <- FALSE # if maxIterations > 1 set allowGroupRemoval <- TRUE


#
# Load sequence from reference data-base: ----
#
# Read sequences into memory
seqs <- readDNAStringSet(seqs_path) # or  readRNAStringSet for RNA
# if gaps in sequences:
seqs <- RemoveGaps(seqs)

#
# Read the taxonomic assigments from data-base ----
# 
taxonomy <- read.table(tax_path, header = FALSE, stringsAsFactors = FALSE)

# Renombramos el nombre de secuencias con la asignacion correspondiente:
if (identical(taxonomy$V1, names(seqs))) {
               names(seqs) <- taxonomy$V2
               groups <- names(seqs)
               groups <- gsub("(.*)(Root;)", "\\2", groups)
               groupCounts <- table(groups)
               cat("number of groups in reference are:\n", length(u_groups <- names(groupCounts)))
    } else {
               groups <- taxonomy$V2
               groups <- gsub("(.*)(Root;)", "\\2", groups)
               groupCounts <- table(groups)
               cat("number of groups in reference are:\n", length(u_groups <- names(groupCounts)))
               }
                 
#
cat('n/', 'Process the training set of data-base', '\n')
#

# Pruning the traning set
# Count the number of representative per group

# maxGroupSize <- 2 # max sequences per label (>= 1); can be set to Inf (infinity) to allow for an unlimited number of sequences per group.
remove <- logical(length(seqs))

for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index),
                 maxGroupSize)
  remove[index[-keep]] <- TRUE
}

save.image(file = paste0(path,'/', out_name, '.RData'))

cat('\n','number of sequences eliminated', sum(remove),'\n') 


```

Debido a que la tarea es computacionalmente costosa, evaluamos la siguiente etapa en otro script:

```R
# Iteratively training the classifier
# Is it will identify any training sequences whose 
# assigned classifications completely (with very high confidence) disagree with their predicted classification.

# maxIterations <- 1 # must be >= 1; also can remove mislabeled in the training data by turning itinerations > 1 
# allowGroupRemoval <- FALSE

library(DECIPHER)

probSeqsPrev <- integer()
taxid <- NULL

load('run015-Mock27-COI-Zoo_OMEGA_A_1e_120_asvs_BOLD_public_species_DECIPHER.taxonomy.RData')

for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove],
                           names(seqs)[!remove],
                           taxid)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}


# saving learn-taxa step: ----
save(trainingSet,
     out_name,
     path,
     q_path,
     file = paste0(path,'/', out_name, '.trainingSet.RData'))

# OR
saveRDS(trainingSet, file = "trainingSet.rds")

quit(save = 'no')
```

Finalmente usamos el objeto `trainingSet` para clasificar nuestras secuencias de interes: 

```R

# View training results:

#trainingSet
#plot(trainingSet)

#
# Classifying Sequences: ----
# 

library(DECIPHER)

f <- 'run015-Mock27-COI-Zoo_OMEGA_A_1e_120_asvs_BOLD_public_species_DECIPHER.taxonomy.trainingSet.RData'

newenv <- new.env()

load(file=f, env=newenv)

q_path <- newenv$q_path

# this is the gold-object for this stem
trainingSet <- newenv$trainingSet 

                           
out_name <- newenv$out_name
path <- newenv$path

processors <- 1
threshold <- 0 # Each taxonomic level is given a confidence between 0% and 100%

q_seqs <- readDNAStringSet(q_path)
q_seqs <- RemoveGaps(q_seqs) # if gaps remove it

ids <- IdTaxa(q_seqs, 
              trainingSet, 
              type="extended",
              strand = "both",
              threshold = threshold,
              processors = processors)

#ids[1:5]
#ids[[2]]
#ids[c(10, 25)]
#c(ids[10], ids[25])

plot(ids, trainingSet)

# subset any rank
rank <- sapply(ids,
                 function(x) {
                   w <- which(x$rank=="Root")
                   if (length(w) != 1) {
                     "unknown"
                   } else {
                     x$taxon[w]
                   }
                 })
table(rank)

taxon <- sapply(ids,
                  function(x)
                    x$taxon[length(x$taxon)])

# Output assignments: ----
output <- data.frame(sapply(ids,
                     function(x)
                       paste(x$taxon,
                             collapse=";")))


output <- data.frame(ASV = rownames(output), Taxonomy = output[,1])
                            
write.table(output, paste0(path, '/', out_name,'2'), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      

 strsplit(output_confidence[, 9], ";")                      
# boots <- as.data.frame(apply(boots0, 2, function(x) gsub("[A-z||()]", "",  x, perl=TRUE)), stringsAsFactors = F)
                            
output_confidence <- sapply(ids,
                 function (id) {
                   paste(id$taxon,
                         " (",
                         round(id$confidence, digits=1),
                         sep="",
                         collapse="; ")})
                            
                            
output_confidence <- as.data.frame(output_confidence, stringsAsFactors = FALSE)
                            
conf_split <- strsplit(output_confidence[,], ";") 
max.rank <- max(lengths(conf_split))
boots0 <- sapply(conf_split, "[", c(1:max.rank))
boots0 <- as.data.frame(t(boots0))

boots <- as.data.frame(apply(boots0, 2, function(x) gsub("[A-z||()]", "",  x, perl=TRUE)), stringsAsFactors = F)
                 
boots[is.na(boots)] <- 100
boots <- apply(boots, 2, as.numeric) # check NAs
boots <- data.frame(boots)


writeLines(output_confidence, paste0(path, '/', out_name, ".confidence"))


quit(save='no')

```



## Testing training set with multiple parameters

...

```scp rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/PRUEBAS/IDTAXA_decipher/trainingSet.rds .```

```R
library(DECIPHER)

#
source(file = "~/Documents/GitHub/metagenomics/readtx.R")

q_tax <- 'BOLD_public.ALL.test.tax'

taxonomy.obj <- read.csv(q_tax, header=FALSE, sep="\t", stringsAsFactors=FALSE)

tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")

max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
taxonomy <- as.data.frame(t(taxonomy))

# Restore the object

f <- 'run015-Mock27-COI-Zoo_OMEGA_A_1e_120_asvs_BOLD_public_species_DECIPHER.taxonomy.trainingSet.RData'

newenv <- new.env()

load(file=f, env=newenv)


## Special case
#path <- '/Users/cigom/metagenomics/db/bold/BOLD_public_trim'
q_path <- 'BOLD_public.ALL.test.fasta'


# this is the gold-object for this stem
trainingSet <- newenv$trainingSet 

                           
out_name <- newenv$out_name
path <- newenv$path

processors <- 3
threshold <- 0 





q_seqs <- readDNAStringSet(q_path)
q_seqs <- RemoveGaps(q_seqs) # if gaps remove it

processors <- 2
threshold <- 100 # Each taxonomic level is given a confidence between 0% and 100%

ids <- IdTaxa(q_seqs, 
              trainingSet, 
              type="extended",
              strand = "both",
              threshold = threshold,
              processors = processors)

taxon <- sapply(ids,
                  function(x)
                    x$taxon[length(x$taxon)])
               
v2 <- data.frame(table(taxon))
                
v1 <- data.frame(table(taxonomy$V9))

                
v1[v1$Var1 %in% v2$taxon, ]

```



### Concepts

> Ref: https://drive5.com/usearch/manual8.1/tax_err.html

**False positive error**: Classifier predicts an incorrect taxon for the given level (family, genus etc.).

**Misclassification error **Type of false positive error. Classifier predicts a wrong name at the given level (family, genus etc.) when at least one reference sequence for the correct taxon is present in the training set.

**Overclassification error** Type of false positive error. Classifier predicts a name at the given level (family, genus etc.) when there are no reference sequences for the correct taxon in the training set. See [taxonomy overclassification and underclassification errors](https://drive5.com/usearch/manual8.1/tax_overclass.html) for further discussion.

**False negative error** Classifier does not predict a name for the given taxon (family, genus etc.) when at least one reference sequence for the correct taxon is in the training set.

**Underclassification error **Classifier does not predict a name at the given level (family, genus etc.) when there are reference sequences for the correct taxon in the training set. See [taxonomy overclassification and underclassification errors](https://drive5.com/usearch/manual8.1/tax_overclass.html) for further discussion. All false negatives are underclassification errors so there is really no need for a new term.

| **Taxon is present in training data** | **Taxon is predicted by classifier** | **Correct name is predicted** | **Correct / Error** |              **Result**              |
| ------------------------------------- | ------------------------------------ | ----------------------------- | ------------------- | :----------------------------------: |
| Yes                                   | Yes                                  | Yes                           | Correct             |            True positive             |
| Yes                                   | Yes                                  | No                            | Error               |  False positive (misclassification)  |
| Yes                                   | No                                   | *-*                           | Error               | False negative (underclassification) |
| No                                    | Yes                                  | *-*                           | Error               | False positive (overclassification)  |
| No                                    | No                                   | *-*                           | Correct             |            True negative             |

 