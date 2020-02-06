3. **RDPtools** ([@Ricardo Gore](https://app.asana.com/0/455610172835450/list))

*- Definir confiabilidad asignacion / CURVAS ROC de la clasificacion (rdp) y procesamiento en R caret (usar subset de datos p68 extraido de bold/gb)*

*- que requiere rdp de input? como procesar los output para curvas roc en r?* 

# Training the classifier

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

`classifier taxa-sim samplefiles/new_trainset_db_taxid.txt samplefiles/new_trainset.fasta samplefiles/Armatimonadetes.fasta taxa_sim 8 rankFile.txt sab`

```bash
#!/bin/bash
#SBATCH --job-name=taxSim
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

rdp=/LUSTRE/bioinformatica_data/genomica_funcional/bin/rdptools/share/rdptools-2.0.2-1/

export PATH=$rdp:$PATH

queryFile=$1
trainSeqFile=$2
trainTaxonFile=$3
rankFile=$4

out=`basename ${queryFile%.fasta}`

mkdir -p ${out}_taxa_sim

classifier taxa-sim $trainTaxonFile $trainSeqFile $queryFile ${out}_taxa_sim 8 $rankFile sab

exit

# taxonfile trainset.fasta query.fasta outdir kmersize rankFile sab|pw
```

How to run in the cluster

```bash
sbatch taxa_sim.sh 
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



## Convert formats 

