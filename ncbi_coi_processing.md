# NCBI Entrez Direct UNIX E-utilities 

It is available to download from the NCBI website [here](ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect) or [here](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

Dowload full ncbi COI barcodes data-base using ncbi nucleotide search as Teresita et.al (2018):

>  CO1[GENE] OR COI[GENE] OR COX1[GENE] OR COXI[GENE] AND "Eukaryota"[Organism] 



- Download Date: July/9/2019

- Number of sequence: 3,388,817 (2,099,002 uniques)
- Size : **5.2G** 
- Filename: CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta
- **Cluster path**: /LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/ncbi/CO1_COI_COX1_COXI_GENE_Eukaryota_nr/CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta
- Downolad method: **e-utilities**
- 

```bash
esearch -db nucleotide -query "CO1[GENE] OR COI[GENE] OR COX1[GENE] OR COXI[GENE] AND "Eukaryota"[Organism] " | efetch -format fasta > CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta
```

Further details of COI barcodes

```bash
# XML format
esearch -db nucleotide -query "CO1[GENE] OR COI[GENE] OR COX1[GENE] OR COXI[GENE] AND "Eukaryota"[Organism] " | efetch -db  taxonomy -format xml > CO1_COI_COX1_COXI_GENE_Eukaryota_nr.xml

# get lineage 
grep 'OrgName_lineage' CO1_COI_COX1_COXI_GENE_Eukaryota_nr.taxonomy

# get 
grep 'OrgMod_subname'  CO1_COI_COX1_COXI_GENE_Eukaryota_nr.taxonomy
```



>  Additional steps [here](http://bioinformatics.cvr.ac.uk/blog/ncbi-entrez-direct-unix-e-utilities/)



## Taxonomy

Get the lineage of the sequences:

1) 

```bash
# 1 get id
grep '^>' CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta  | awk '{print $1}' | sed 's/>//' > ids

# 2 search id lineage
for i in $(cat ids); do esearch -db nucleotide -query $i | efetch -format xml | xtract -pattern Org-ref -element Object-id_id,Org-ref_taxname,OrgMod_subname,OrgName_lineage,Textseq-id_accession,OrgName_div,BinomialOrgName_genus,BinomialOrgName_species; done > CO1_COI_COX1_COXI_GENE_Eukaryota_nr.xtract.Org-ref-element
```

2) 

```bash
# full lineage w/ taxid
esearch -db protein -query "NP_066243"| elink -target taxonomy |efetch -format xml |xtract -pattern TaxaSet -element Lineage,Rank,TaxId
```



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

> primer GGWACWGGWTGAACWGTWTAYCCYCC TAIACYTCIGGRTGICCRAARAAYCA

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