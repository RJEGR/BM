# Transcriptome Assembly Quality Assessment

The aim of (denovo) transcriptome assembly is to accurately reconstruct the complete set of transcripts that are represented in the read data (in the absence of reference genome)

In contrast of genome assembly, (**context**: if we have an organism with C chromosomes our optimal (ideally) genome assembly would consist of C long contigs) transcrimptome assembly (optimal) wil vary consist of all posible transcript from all expressed genes. ie. alternatively spliced variants (isoforms). 

Nevertheless there are several contributing factors than negatively affect the accuracy of the construction assembly process.

These factors include error in sequencing process (usually the batch effect is found here), incomplete coverage of transcripts (due to insuficient sequencing depth), real biological variability (ex. alternatively spliced variants) and algorithmic simplification. 

Usualy the most highly expressed transcript do not neccessary constitute the longjest one and the majority of transcripts in a transcriptome assembly will normally have relatively low expression levels.



# N50

Imagine that you line up all the contigs in your assembly in the order of their sequence lengths (Fig. 1a). You have the longest contig first, then the second longest, and so on with the shortest ones in the end. Then you start adding up the lengths of all contigs from the beginning, so you take the longest contig + the second longest + the third longest and so on — all the way until you’ve reached the number that is making up 50% of your total assembly length. That length of the contig that you stopped counting at, this will be your N50 number (Fig. 1b).

We strongly advise against using regular N50 metrics for transcriptome assemblies. Instead, other more appropriate measures can be used. The developers of the transcriptome assembler Trinity have invented [the **ExN50** metric](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats), which takes into account the expression levels of each contig and is therefore a more suitable contig length metric for transcriptomes.  

# Ex90

...

# Reads represented

To assess the read composition of our assembly, we want to capture and count all reads that map to our assembled transcripts, including the properly paired and those that are not we run the process below. Bowtie2 is used to align the reads to the transcriptome and then we count the number of proper pairs and improper or orphan read alignments. First, build a bowtie2 index for the transcriptome and then perform the aligment to catpure the paired-reads aligment statistic. (Ref [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly))



```bash


bowtie2-build Trinity.fasta Trinity.fasta
# Then perform the alignment (example for paired-end reads) to just capture the read alignment statistics.

srun bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 R1.P.qtrim.fq -2 R2.P.qtrim.fq | samtools view -@10 -Sb -o ./bowtie2.bam


```

Finally lets summary the bam file using bamtools:

```bash
bamtools stats -in reads_represented/bowtie2.bam

#or
util=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/

export PATH=$util:$PATH

SAM_nameSorted_to_uniq_count_stats.pl reads_represented/bowtie2.bam

# or
support_scripts=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/support_scripts/

export PATH=$support_scripts:$PATH

srun SAM_coordSorted_fragment_Read_coverage_writer.pl bowtie2.bam &> SAM_coordSorted_fragment_Read_coverage_writer.log &

# 
srun SAM_extract_properly_mapped_pairs.pl bowtie2.bam &> SAM_extract_properly_mapped_pairs.txt &

# and check SAM_ scriptss
SAM_extract_uniquely_mapped_reads.pl

srun SAM_extract_uniquely_mapped_reads.pl bowtie2.bam &> SAM_extract_uniquely_mapped_reads.txt &

# Get ids
# 142,725 ids
grep 'TRINITY_' SAM_extract_uniquely_mapped_reads.txt | awk '{print $3}' | sort | uniq > SAM_extract_uniquely_mapped_reads.ids
# 361,343 ids
grep 'TRINITY_' SAM_extract_properly_mapped_pairs.txt | awk '{print $3}' | sort | uniq > SAM_extract_properly_mapped_pairs.ids

# Get sequences 
cat Trinity.fasta | seqkit grep -f SAM_extract_uniquely_mapped_reads.ids > Trinity_SAM_uniquely_mapped_reads.fasta

cat Trinity.fasta | seqkit grep -f SAM_extract_properly_mapped_pairs.ids > Trinity_SAM_extract_properly_mapped_pairs.fasta

# Get stats
TrinityStats.pl Trinity_SAM_uniquely_mapped_reads.fasta
TrinityStats.pl Trinity_SAM_extract_properly_mapped_pairs.fasta


```

> in case you neet, fist install bamtools as follow: https://github.com/pezmaster31/bamtools/wiki/Building-and-installing

A output in your screen will be printed as follow:

**********************************************
Stats for BAM file(s):
**********************************************

Total reads:       13970909
Mapped reads:      13970909	(100%)
Forward strand:    6885283	(49.283%)
Reverse strand:    7085626	(50.717%)
Failed QC:         0	(0%)
Duplicates:        0	(0%)
Paired-end reads:  13970909	(100%)
'Proper-pairs':    11803698	(84.4877%)
Both pairs mapped: 13127153	(93.9606%)
Read 1:            7037955
Read 2:            6932954
Singletons:        843756	(6.03938%)

**********************************************

The [Integrative Genomics Viewer](http://software.broadinstitute.org/software/igv/) is useful for visualizing read support across any of the Trinity assemblies. The bowtie2 alignments generated above, which are currently sorted by read name, can be re-sorted according to coordinate, indexed, and then viewed along with the Trinity assemblies using the IGV browser as follows.

```bash
# sort the alignments by coordinate
samtools sort bowtie2.bam -o bowtie2.coordSorted.bam

# index the coordinate-sorted bam file
samtools index bowtie2.coordSorted.bam

# index the Trinity.fasta file
samtools faidx Trinity.fasta

# view the aligned reads along the Trinity assembly reference contigs.
# note, you can do this by using the various graphical menu options in IGV (load genome 'Trinity.fasta', load file 'bowtie2.coordSorted.bam'), or you can use the command-line tool like so:

igv.sh -g `pwd`/Trinity.fasta  `pwd`/bowtie2.coordSorted.bam
```

## Completness

> example in `/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/oyster_gon/assembly_check/completness`

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=BUSCOpy
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

module load python-2.7-anaconda

fasta=$1
out=${fasta%.fasta}

BUSCO=/home/rgomez/bin/busco-master/scripts/
export PATH=$BUSCO:$PATH

eukaryota_odb9=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/eukaryota_odb9
metazoa_odb10=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/metazoa_odb10
mollusca_odb10=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/mollusca_odb10
bacteria_odb9=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/bacteria_odb9


run_BUSCO.py -i $fasta -l $eukaryota_odb9 -m transcriptome -o ${out}_eukaryota_odb9 -c 24
run_BUSCO.py -i $fasta -l $metazoa_odb10 -m transcriptome -o ${out}_metazoa_odb10 -c 24
run_BUSCO.py -i $fasta -l $mollusca_odb10 -m transcriptome -o ${out}_mollusca_odb10 -c 24
run_BUSCO.py -i $fasta -l $bacteria_odb9 -m transcriptome -o ${out}_bacteria_odb9 -c 24

exit
```

Entonces generamos las figuras

```bash
mkdir summaries
cp ./run_Trinity_*_odb*/short_summary* summaries

module load R-3.3.1

python2.7 /home/rgomez/bin/busco-master/scripts/generate_plot.py --working_directory ./summaries/
cp summaries/busco_figure.R .
sed -i 's/_odb9//g' summaries/busco_figure.R
Rscript summaries/busco_figure.R
firefox summaries/busco_figure.png
cp summaries/busco_figure.R .
```

Tambien retenemos contigs en los que su status de identidad resulto `completo`, `duplicado` y `fragmentado`.

```bash
# 895197 n contigs in full- assembly oktopus

egrep 'Complete|Duplicated|Fragmented' full_table_*.tsv | awk '{print $3}' | sort | uniq > SC_SD_F_TRINITY_IDS.txt

# 2.
# linearizeFasta.awk
# https://bioinf.shenwei.me/seqkit/download/

cat Trinity.fasta | seqkit grep -f run_Trinity_metazoa_odb10/SC_SD_F_TRINITY_IDS.txt > Trinity_metazoa_odb10.fasta

# Metazoa
cat Trinity.fasta | seqkit grep -f run_Trinity_metazoa_odb10/SC_SD_F_TRINITY_IDS.txt > Trinity_metazoa_odb10.fasta

# Mollusca
cat Trinity.fasta | seqkit grep -f run_Trinity_mollusca_odb10/SC_SD_F_TRINITY_IDS.txt > Trinity_mollusca_odb10.fasta

# Eukaryota
cat Trinity.fasta | seqkit grep -f run_Trinity_eukaryota_odb9/SC_SD_F_TRINITY_IDS.txt > Trinity_eukaryota_odb9.fasta

TrinityStats.pl Trinity.fasta
TrinityStats.pl Trinity_metazoa_odb10.fasta
TrinityStats.pl Trinity_eukaryota_odb9.fasta
TrinityStats.pl Trinity_mollusca_odb10.fasta
```

## Counting Numbers of Expressed Transcripts or Genes

```bash
#ls ./LOF_*/RSEM.isoforms.results > isoform.results
 
misc=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/misc
# by gene level

$misc/count_matrix_features_given_MIN_TPM_threshold.pl RSEM.trans.gene.TPM.not_cross_norm | tee genes_matrix.TPM.not_cross_norm.counts_by_min_TPM
     
# gene level ----
assembly=Trinity.fasta
# should down cause gene level ids
$misc/contig_ExN50_statistic.pl RSEM.trans.gene.TMM.EXPR.matrix $assembly | tee genes_ExN50.stats


# and by transcript level ----

$misc/count_matrix_features_given_MIN_TPM_threshold.pl RSEM.trans.isoform.TPM.not_cross_norm | tee trans_matrix.TPM.not_cross_norm.counts_by_min_TPM

assembly=Trinity.fasta

$misc/contig_ExN50_statistic.pl RSEM.trans.isoform.TMM.EXPR.matrix $assembly | tee trans_ExN50.stats
```

Also run script by:

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=QCA
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# Quality Check Assembly methods by trinity: 

misc=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/misc

TPM=$1 # RSEM.trans.gene.TMM.EXPR.matrix
assembly=$2 # Trinity.fasta
cross_norm=$3 # RSEM.trans.gene.TPM.not_cross_norm
# by gene level

$misc/count_matrix_features_given_MIN_TPM_threshold.pl $cross_norm | tee ${cross_norm%.not_cross_norm}_by_min_TPM

$misc/contig_ExN50_statistic.pl $TPM $assembly | tee ${TPM%.matrix}_ExN50.stats

exit

```



### Longest and Abundant isoform

Dejemos de trabajar con isoformas y usemos genes, hay varias formas de trabajar:

- Usar la isoforma mas abundante (`$TRINITY_HOME/util/filter_low_expr_transcripts.pl`) como la representativa

```bash
assembly=Trinity.fasta

sbatch abundance.sh -t $assembly -s sampleFile

#  -------
#!/bin/sh
## Directivas
#SBATCH --job-name=filAb
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# ls ./*/*isoforms.results > isoforms.results

assembly=$1 # Trinity.fasta
results=$2 # isoforms.results
# 1)
# Get the gene-trans map file
support_scripts=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/support_scripts/

$support_scripts/get_Trinity_gene_to_trans_map.pl $assembly > ${assembly%.fasta}.gene_trans_map

# 2) Convert abundance results a matrix
util=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/

$util/abundance_estimates_to_matrix.pl \
        --gene_trans_map ${assembly%.fasta}.gene_trans_map \
        --est_method RSEM \
        --out_prefix RSEM \
        --quant_files $results \
        --name_sample_by_basedir

# 3) Filter low exp transcripts
# this would ideally be your TPM matrix - or TMM-normalized TPM matrix *not* raw counts)

$util/filter_low_expr_transcripts.pl --matrix RSEM.isoform.TMM.EXPR.matrix --transcripts $assembly --min_expr_any 1 --gene_to_trans_map ${assembly%.fasta}.gene_trans_map > ${assembly%.fasta}.min1TPM.fasta

TrinityStats.pl ${assembly%.fasta}.min1TPM.fasta > ${assembly%.fasta}.min1TPM.stats

exit

# continue with stasts: in 
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#filtering-transcripts
```



- Usar la isoforma mas larga como la representativa `$TRINITY_HOME/util/misc/get_longest_isoform_seq_per_trinity_gene.pl`
- Elaborar un scafold de todas las isoformas (SuperTranscript)

```bash
find "$(pwd)" -name '*R2*P.qtrim.gz' | sort > R2_P.qtrim.list
find "$(pwd)" -name '*R1*P.qtrim.gz' | sort > R1_P.qtrim.list


find -name '*R2*P.qtrim.gz' | sort | cut -d'_' -f1 | sed 's/[./]//g' > factor1.list
find -name '*R2*P.qtrim.gz' | sort | cut -d'_' -f1,2,3 | sed 's/[./]//g' > factor2.list

paste -d'\t' factor2.list factor2.list R1_P.qtrim.list R2_P.qtrim.list > samples.file

```



## 

Check /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oktopus_full_assembly/reads_represented/BUSCO

```bash
ls ./LOF_*/RSEM.isoforms.results > isoform.results
 
util=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/

$util/abundance_estimates_to_matrix.pl \
        --gene_trans_map trinity_mapped_pairs.fasta.gene_trans_map \
        --est_method RSEM \
        --out_prefix iso \
        --quant_files isoform.results \
        --name_sample_by_basedir
        

# sed -i 's/_R1//g' iso.counts.matrix

TOOL=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/

module load R-3.5.0_bio
# que no este cargado el source de gmelendrez o no correra
$TOOL/run_DE_analysis.pl \
      	  --matrix iso.counts.matrix \
       	 --method DESeq2 \
       	 --samples_file LOF_samples.file \
       	 --output DESeq2_dir
       	 
       	 #--contrasts contrast_file.txt \

# Filter DiffExp up/down
DF=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression/

$DF/analyze_diff_expr.pl --matrix ../iso.TMM.EXPR.matrix --samples ../LOF_samples.file -P 1e-3 -C 1 --order_columns_by_samples_file

wc -l diffExpr.P*.matrix
wc -l *UP.subset

# 
mkdir stast
cd stast

sbatch PtR.sh iso.counts.matrix LOF_samples.file
# /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oktopus_full_assembly/DiffExp/Diana_results/stats


#



```



Continue with phylogenomics

https://gitlab.com/ezlab/busco_usecases/-/tree/master/phylogenomics

