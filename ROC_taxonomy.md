3. **RDPtools** ([@Ricardo Gore](https://app.asana.com/0/455610172835450/list))

- Probar rdp tools con set de muestra (Hecho)
  - taxa-sim
  - Leave-and-out test
- Dar formato a archivos con datos personales (P68 + gb subset) (Ahora)
  - taxonfile 
  - trainset.fasta 
  - query.fasta
- Procesar datavis en R (EN PROCESO)
  - caret
  - Ggplot2 (HECHO)

# 1. Training the classifier

Follow these steps when there is a need to retrain Classifier, such as novel lineages, newly named type organisms, taxonomic rearrangements, better training set covering specific taxa, or alternative taxonomy. Two files, a taxonomy file and a training sequence file with lineage are required. Prefer high quality, full length sequences, or at least covering the entire region of gene of interest. See samplefiles for example data files.  Based on our experience, trimming the sequences to a specific region does not improve accuracy. The ranks are not required to be uniform neither, which means you can define any number of ranks as necessary. The speed of the Classifier is proportional to the number of genera, not the number of training sequences.

**1. Plot intra taxon Similarity by fraction of matching 8-mer**

Use subcommand "taxa-sim" to calculate and plot intra taxon Similarity by fraction of matching 8-mer (see example plots using fungal ITS training sets on RDP's poster http://rdp.cme.msu.edu/download/posters/MSA2014_RDP.pdf). To run taxa-sim in Headless mode without GUI display, use the following options:

> `rankFile`: a file contains a list of ranks to be calculated and plotted. One rank per line, no particular order required:
>
> ```
> domain
> phylum
> class
> order
> family
> genus
> species
> ```

```bash
#!/bin/bash
#SBATCH --job-name=taxSim
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

rdp=/LUSTRE/bioinformatica_data/genomica_funcional/bin/rdptools/share/rdptools-2.0.2-1/

queryFile=$1
trainSeqFile=$2
trainTaxonFile=$3
rankFile=$4

out=`basename ${queryFile%.fasta}`

mkdir -p ${out}_taxa_sim

java -Djava.awt.headless=true -jar $rdp/classifier.jar taxa-sim $trainTaxonFile $trainSeqFile $queryFile ${out}_taxa_sim 8 $rankFile sab

exit

# taxonfile trainset.fasta query.fasta outdir kmersize rankFile sab|pw
```

How to run in the cluster:

```bash
sbatch taxa_sim.sh samplefiles/Armatimonadetes.fasta samplefiles/new_trainset.fasta samplefiles/new_trainset_db_taxid.txt rankFile.txt
```



**2. Estimate the accuracy of your own training data using leave-one-out testing:**

The program will output a tab-delimited test result file which can be loaded to Excel and plot the accuracy rates. It also contains the list of misclassified sequences and the rank when misclassified seqs group by taxon. Examine the result careful to spot errors in the taxonomy.

> <u>Inputs</u> 
>
> -q --queryFile query file contains sequences, same format as the training sequence file
> -s --trainSeqFile training files in fasta format labelled with the lineage information
>
> The header of this fasta file starts with '>', followed by the sequence name, white space(s) and a list taxon names seperated by ';' with highest rank taxon first ex: Root;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter

> ```
> >EF515962|S000840404  Root;Bacteria;"Armatimonadetes";Chthonomonadetes;Chthonomonadales;Chthonomonadaceae;Chthonomonas/Armatimonadetes_gp3
> acgaacgcttgcggcgtgcctaagaaatgcaagtcgagcggggagcaat ...
> ```
>
> -t --trainTaxonFile contains the hierarchical taxonomy information, taxon name and rank together is unique. The format looks like the following: taxid*taxon name*parent taxid*depth*rank Note taxid, the parent taxid and depth should be in integer format. depth indicates the depth from the root taxon. Recommend removing duplicate seqeunces using command `rmdupseq`.
>
> ```
> 0*Root*-1*0*rootrank
> 1*Bacteria*0*1*domain
> 2*"Actinobacteria"*1*2*phylum
> 3*Actinobacteria*2*3*class
> 4*Acidimicrobidae*3*4*subclass
> 5*Acidimicrobiales*4*5*order
> 6*"Acidimicrobineae"*5*6*suborder
> 7*Acidimicrobiaceae*6*7*family
> 8*Acidimicrobium*7*8*genus
> 9*Ferrimicrobium*7*8*genus
> 10*Ferrithrix*7*8*genus
> ```

**a) Leave-one-sequence-out testing:** 

Each iteration one sequence from the training set was chosen as a test sequence. That sequence was removed from training set. The assignment of the sequence produced by the Classifier was compared to the original taxonomy label to measure the accuracy of the Classifier.

```bash
#!/bin/bash
#SBATCH --job-name=loot
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

rdp=/LUSTRE/bioinformatica_data/genomica_funcional/bin/rdptools/share/rdptools-2.0.2-1/

export PATH=$rdp:$PATH

queryFile=$1
trainSeqFile=$2
trainTaxonFile=$3

out=`basename ${queryFile%.fasta}`

# loso  - Leave-one-sequence-out testing: 
classifier loot -q $queryFile -s $trainSeqFile -t $trainTaxonFile -l 400 -o ${out}_loso.txt

# loto - Leave-one-taxon-out testing: 
# -h --hideTaxon If set, remove the lowest taxon where a query sequence originally labelled from the training set. Default only remove the query seq from training set

classifier loot -h -q $queryFile -s $trainSeqFile -t $trainTaxonFile -l 400 -o ${out}_loto.txt

exit
```

**b) Leave-one-taxon-out testing:** 

Similar to the leave-one-sequence-out testing except for each test sequence, the lowest taxon that sequence assigned to (either species or genus node) was removed from the training set. This is intended to test if the species or genus is no present in the training set, how likely the Classifier can assign the sequence to the correct genus or higher taxa.

```bash
#!/bin/bash
#SBATCH --job-name=loto
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

rdp=/LUSTRE/bioinformatica_data/genomica_funcional/bin/rdptools/share/rdptools-2.0.2-1/

export PATH=$rdp:$PATH

queryFile=$1
trainSeqFile=$2
trainTaxonFile=$3

out=`basename ${queryFile%.fasta}`

classifier loot -h -q $queryFile -s $trainSeqFile -t $trainTaxonFile -l 400 -o ${out}_loto.txt

exit
```

To run in the cluster

```bash
sbatch rdp_loot.sh samplefiles/Armatimonadetes.fasta samplefiles/new_trainset.fasta samplefiles/new_trainset_db_taxid.txt
```

Results available are:

```bash
# grep '^[**]' Armatimonadetes_loso.txt
**The statistics for each hierarchy level:
** 2. The average votes for each bin range
** 3. The percentage of correctness for each bin range (the percentage of #1)
** 4. The standard error for each bin range
**misclassified sequences:
**singleton sequences:
**misclassified sequences group by taxon
**ROC matrix
**Area under curve
```

Split file by repex

```bash
csplit -z Armatimonadetes_loso.txt '/^[**]/' '{*}' # in bash
mv xx08 ROC_Armatimonadetes_loso.txt
rm xx*
#
csplit -z Armatimonadetes_loto.txt '/^[**]/' '{*}' # in bash
mv xx08 ROC_Armatimonadetes_loto.txt
rm xx*
```



## 2. Process results in R

ROC matrix

```

s <- read.table(file1, sep="\t", skip = 2, header = T)
t <- read.table(file2, sep="\t", skip = 2, header = T)

bootstrap	
rank_FPR	
rank_TPR	
rank_F1score
...		
genus_FPR	
genus_TPR	
genus_F1score	

```



## Convert formats

> MD de M. Martinez, 2020

1. Sequence file (e.g. rawSeqs.fasta)

  * in fasta format with unique identifier (string without spaces) for each sequence 

2. Taxonomy file (e.g. rawTax.txt)
  * txt, TAB-sep

  ```bash
  cat ncbi_complete_drep_subset.taxonomy | awk '{gsub(";","\t"); print $1,"Root", $2, $3, $4, $5, $6, $7}'| column -t > ncbi_complete_subset.taxonomy
  ```

  * Header:  First column: sequence identifier; the following columns contain taxonomic rank names one in a column in the order from root (highest) to leaf rank (lowest), such as Domain/Kingdom, Phylum, Class,Order, Family, Genus, etc.) for each  taxon level you want to represent:

  ```bash
  vi ncbi_complete_drep_subset.convergent.taxonomy
  # acc	Kingdom	Class	Order	Family	Genus	Specie
  # id	Root	Domain	Kingdom	Phylum	Class	Order	Family
  ```

  * Columns should be:
  	1. SeqID.  Same as in fasta file
  	2-N Rank Names. in order from higher (kingdom) to lower (genus)
  	
  	```bash
  	0. AB477020.1
  	1. root
  	2. Metazoa
  	3. Chordata  
  	4. Actinopteri  
  	5. Pleuronectiformes
  	6. Cynoglossidae     
  	7. Paraplagusia
  	```
  	
  * Empty taxon names shuold be '-'

  * >  **Convergent taxons are not allowed** (same taxon name for two different parents, e.g. Clostridiaceae>Clostridium and Eubacteria>Clostridium

  ```bash
  file=ictio_coi_sanger_bold_trust.fasta # ncbi_complete_drep.fasta
  grep '^>' $file | sed 's/>//g' | head -n100 > ${file%.fasta}_subset.id
  
  # Get taxonomy based on a list File
  
  tax=ictio_coi_sanger_bold_trust.tax # ncbi_complete.taxonomy
  
  rm ${tax%.tax}_subset.tax
  
  while IFS= read -r pattern; do grep $pattern $tax | awk '{gsub(";" , "\t"); print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6, "\t",$7,"\t",$8}' >> ${tax%.tax}_subset.tax  & done < ${file%.fasta}_subset.id; wait; echo "Grouping  genus done"
  
  # Get unique taxons to clean convergence
  # awk '{print $8}' ${tax%.tax}_subset.tax | sort
  
  awk '{print $8}' ${tax%.tax}_subset.tax | sort | uniq > ${tax%.tax}_subset.convergent
  
  # Clean convergent taxons
  
  rm ncbi_complete_drep_subset.convergent.taxonomy
  
  while IFS= read -r pattern; do grep -w $pattern ncbi_complete_drep_subset.taxonomy | head -n1 >> ncbi_complete_drep_subset.convergent.taxonomy  & done < ncbi_complete_drep_subset.convergent; wait; echo "Grouping sp"
  
  #awk -F":" -v OFS=',' ncbi_complete_drep_subset.convergent
  # awk 'BEGIN {FS="\t"}; {print $0}' ncbi_complete_drep_subset.convergent
  
  wc -l ncbi_complete_drep_subset.convergent.taxonomy
  
  # sanity check 
  
  awk '{print $8}' ncbi_complete_drep_subset.convergent.taxonomy | sort | uniq | wc -l
  
  # 
  cut -f1 ncbi_complete_drep_subset.convergent.taxonomy > ncbi_complete_drep_subset.convergent.id
  
  # Get fasta (not working!!!!)
  cat ncbi_complete_drep.fasta | seqkit grep -f ncbi_complete_drep_subset.convergent.id > ncbi_complete_drep_subset.convergent.fasta
  
  # (Not working)
  seqkit grep -f ncbi_complete_drep_subset.convergent.id ncbi_complete_drep.fasta > ncbi_complete_drep_subset.convergent.fasta
  ```

  Then run:

  ```bash
  ./lineage2taxTrain.py ncbi_complete_drep_subset.convergent.taxonomy > ncbi_complete_drep_subset.convergent.taxTrain
  ```

  Fields taxid, the parent taxid and depth should be in integer format.  The taxid, or the combination of taxon name and rank is unique depth indicates the depth from the root taxon. Note: the depth for the root is 0

  ```
  0*Root*-1*0*rootrank
  1*Bacteria*0*1*domain
  2*"Actinobacteria"*1*2*phylum
  3*Actinobacteria*2*3*class
  ```

  `lineage2taxTrain.py`

  ```bash
  #!/usr/bin/env python
  #used to convert a taxonomy in tab-delimited file containing the taxonomic hierarchical structure to RDP Classifier taxonomy training file
  
  #Approach:each taxon is uniquely identified by the combination of its tax id and depth from the root rank, its attributes comprise: name, parent taxid, and level of depth from the root rank. 
  
  import sys, string
  
  if not len(sys.argv) == 2:
  	print "lineage2taxTrain.py taxonomyFile"
  	sys.exit()
  
  f = open(sys.argv[1], 'r').readlines()
  header = f[0]
  cols = header.strip().split('\t')[1:]
  hash = {}#taxon name-id map
  ranks = {}#column number-rank map
  lineages = []#list of unique lineages
  
  hash = {"Root":0} #initiate root rank taxon id map
  
  for i in range(len(cols)):
  	name = cols[i]
  	ranks[i] = name
  root = ['0', 'Root', '-1', '0', 'rootrank'] #root rank info
  print string.join(root, '*')
  ID = 0 #taxon id
  
  for line in f[1:]:
  	cols = line.strip().split('\t')[1:]
  #	if not cols in lineages:#unique lineage
  #		lineages.append(cols)
  	for i in range(len(cols)):#iterate each column
  		#name = string.join(cols[:i + 1], ';')
  		name = []
  		for node in cols[:i + 1]:
  			if not node == '-':
  				name.append(node)
  		pName = string.join(name[:-1], ';')
  		if not name in lineages:
  			lineages.append(name)
  		depth = len(name)
  		name = string.join(name, ';')
  		if name in hash.keys():
  			continue
  		rank = ranks[i]
  		#level = len(name.split(';'))
  		#pName = string.join(cols[:i], ';')#parent name
  		if i == 0:
  			pName = 'Root'
  		pID = hash[pName]#parent taxid
  		ID += 1
  		hash[name] = ID #add name-id to the map
  		#out = ['%s'%ID, name, '%s'%pID, '%s'%depth, rank] 
  		out = ['%s'%ID, name.split(';')[-1], '%s'%pID, '%s'%depth, rank] 
  		print string.join(out, '*')
  				
  ```

```bash
classifier merge-detail
# merge-detail  - merge classification detail result files to create a taxon assignment counts file
# merge-count   - merge multiple taxon assignment count files to into one count file

0*Root*-1*0*rootrank
1*Bacteria*0*1*domain
2*"Actinobacteria"*1*2*phylum
3*Actinobacteria*2*3*class
# -----

Root;Bacteria;"Armatimonadetes";Chthonomonadetes;Chthonomonadales;Chthonomonadaceae;Chthonomonas/Armatimonadetes_gp3q

#1. tomar niveles por columnas, y asignar nombres de columnas al rank level respectivo
# 2) por columna asignar nombres unicos 
```

**Prepare fasta file:**

`./addFullLineage.py ncbi_complete_drep_subset.convergent.taxTrain rawSeqs.fasta > ready4train_seqs.fasta`

```bash
#!/usr/bin/python
import sys, string
if len(sys.argv) != 3:
	print 'addFullLineage.py taxonomyFile fastaFile'
	sys.exit()
f1 = open(sys.argv[1], 'r').readlines()
hash = {} #lineage map
for line in f1[1:]:
	cols = line.strip().split('\t')
	lineage = ['Root']
	for node in cols[1:]:
		if not node == '-':
			lineage.append(node)
	ID = cols[0]
	lineage = string.join(lineage, ';')
	hash[ID] = lineage
f2 = open(sys.argv[2], 'r').readlines()
for line in f2:
	if line[0] == '>':
		ID = line.strip().replace('>', '')
		try:
			lineage = hash[ID]
		except KeyError:
			print ID, 'not in taxonomy file'
			sys.exit()
		print line.strip() + '\t' + lineage
	else:
		print line.strip()
```



Vamos a tomar un subconjunto de datos de ncbi.

Estos datos corresponden a 70 ordenes reportados a lo largo de las **89,847** secuencias filtradas en la busqueda en NCBI. Se agrupan ordenes <= 100 secuencias dentro del grupo *Others* (corresponde a 17 / 70 ordenes). Se observa un desfase en 5337 secuencias con orden desconocido "-" pero especie clasificada; esto se debe a los problemas de *gaps* en los linajes de ncbi. En cuanto a la distribución cumulativa de géneros (~ 3169 ). Casi el 99 % de los géneros tienen numeros de lecturas por debajo de los cientos (figura izquierda). Es decir, géneros con secuencias <= 100 (2968 géneros) (figura derecha).



**Vamos a buscar un genero que se nivele en varias especies (cargar datos taxonomicos en R y repasar codigo). Entonces usamos para convertir datos, adicional usamos el set de mock68**

```bash
>EF515962|S000840404  Root;Bacteria;"Armatimonadetes";Chthonomonadetes;Chthonomonadales;Chthonomonadaceae;Chthonomonas/Armatimonadetes_gp3
acgaacgcttgcggcgtgcctaagaaatgcaagtcgagcggggagcaat ...
```

```bash
file = 'ictio_coi_sanger_bold_trust.fasta'
tax <- 'ictio_coi_sanger_bold_trust_sorted.tax'

dna <- readDNAStringSet(file, format="fasta")

x <- read.csv(tax, header=FALSE, sep='\t', stringsAsFactors=FALSE)

x <- apply(x,2,function(x)gsub('\\s+', '',x))
 
tax <- data.frame(x[,2:8])
id <- names(dna)

identical(id, x[,1])

library(tidyr)

save <- cbind(id, unite(tax, sep = ";", remove = TRUE, col = 'Taxonomy'))


save$Taxonomy <- sapply(save$Taxonomy, 
                        function(x){gsub(pattern = "$",
                        replacement = ";", x)})

names <- unite(save, sep = " ", remove = TRUE, col = 'names')

                      
names(dna) <- names$names
                        
writeXStringSet(dna, paste0(getwd(),'/ictio_coi_sanger_bold_trust_header.fasta'), format="fasta")
```

