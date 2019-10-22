Quick-run

```bash
#!/bin/bash
prefix=test
# 1)
sbatch demultiplex.sh raw/ $prefix
# raw directory contain only the raw paired-end files

# 2) after visualize the f and r trimming positions run:
sbatch denoise-paired.sh 150 146 ${prefix}_demux_seqs.qza sampleMfile

# 3) Filter low-frequency features
sbatch filter-features.sh ${prefix}_table.qza 1 sampleMfile

# 4) Classifying features
sbatch Classify-feature.sh ${prefix}_rep-seqs_table_2_filtered.qza 0.8 w2pr2_worms_API02_no_gaps.fasta.qza w2pr2_worms_API02.tax.qza

# 5) Prepare outputs and pre-visualziations
sbatch export-features.sh ${prefix}_table_1_filtered.qza ${prefix}_rep-seqs_table_1_filtered.qza w2pr2_worms_API02_vs_${prefix}_rep-seqs_table_1_filtered.qza

# 6) Write a report in R



```



Explanation script-to-script below.

The framework version of qiime are:

**q2cli version 2018.8.0**

Some modules to pre-load

```bash
module load miniconda3-python-3.5
source activate quiime2-2018.0

# Additional step
conda info --envs
```

## Starting

## Demultiplex

Copy and paste the follow code in a bash script `demultiplex.sh`; then run.

Example: `sbatch demultiplex.sh wordirectory prefix`

The work directory must contain the paired-end fastq files. The prefix is a free-tag to use through the outputs in the analysis with qiime

```bash
#!/bin/bash
#SBATCH --job-name=qiime
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=12
#SBATCH -t 6-00:00:00

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

# Load qiime
module load miniconda3-python-3.5
source activate quiime2-2018.0

# My vars
seqs_dir=$1 # Directory where only libs are found
prefix=$2 # prefix to add in outputs

# my code

# Import reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $seqs_dir \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ${prefix}_demux_seqs.qza

# Visualize results and set trimming parameters

qiime demux summarize \
  --i-data ${prefix}_demux_seqs.qza \
  --o-visualization ${prefix}_demux_seqs.qzv

exit

```

## Denoise



After run `demultiplex.sh` script, lets run the trimming and  denoising step using the follow script

Example: `sbach denoise-paired.sh 150 146 demux_seqs.qza sampleMfile` 

The script `denoise-paired.sh` is as follow

```bash
#!/bin/bash
#SBATCH --job-name=qiime
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# Load qiime
module load miniconda3-python-3.5
source activate quiime2-2018.0

# My vars
trunc_f=$1
trunc_r=$2
demux_seqs=$3
metadt=$4

prefix=$(ls $demux_seqs | cut -d"_" -f1)

qiime dada2 denoise-paired \
	--p-n-threads $SLURM_NPROCS \
  --i-demultiplexed-seqs $demux_seqs \
  --p-trim-left-f 1 \
  --p-trim-left-r 1 \
  --p-trunc-len-f $trunc_f \
  --p-trunc-len-r $trunc_r \
  --o-table ${prefix}_table.qza \
  --o-representative-sequences ${prefix}_rep-seqs.qza \
  --o-denoising-stats ${prefix}_denoising-stats.qza
  
exit

```

## Metadata

The `$metadt` variable must contain first column as SampleID from data. Column of barcode are not requiered fort qiime2. The columns from Factor1 to FactorN are independent . The column separator must be tabular Example:

| #SampleID |      Barcode      | Run  | Cruize | Station | Factor1 | FactorN |
| --------- | :---------------: | :--: | ------ | ------- | ------- | ------- |
| G44.X04   | CTCTCTAC+CGTCTAAT |  9   | X04    | G44     | N       | >4.4    |
| A8.X04    | TACGCTGC+AAGGAGTA |  9   | X04    | A8      | D       | >4.4    |
| H46.X04   | CGAGGCTG+CGTCTAAT |  9   | X04    | H46     | D       | >4.4    |
| Y1A.X04   | GTAGAGGA+CGTCTAAT |  9   | X04    | Y1A     | N       | >4.4    |
| G42.X04   | CGGAGCCT+TCTCTCCG |  9   | X04    | G42     | D       | >4.4    |
| Y4B.X04   | ATGCGCAG+TCTCTCCG |  9   | X04    | Y4B     | D       | >4.4    |
| D26.X04   | ATGCGCAG+CTAAGCCT |  9   | X04    | D26     | N       | >4.4    |
| Y3A.X04   | ATCTCAGG+CGTCTAAT |  9   | X04    | Y3A     | N       | >4.4    |

>  You can use [keemei](https://chrome.google.com/webstore/detail/keemei/ohdhmeoedjeepljniehbmmcgoelhfmop) adds-on to check the metadata-qiime format within google-sheets
>
> Based on the filename you can construct your metadata as follow: 
>
> ```bash
> ls *fastq.gz | cut -d"_" -f1 | uniq | sed 's/-/\t/g' > metadata.tmp
> ls *fastq.gz | cut -d"_" -f1 | uniq > sampleid.tmp
> 
> paste -d"\t" sampleid.tmp metadata.tmp > sampleMfile
> sed -i '1iSampleID\tRun\tCruize\tsample\tmarker\ttype' sampleMfile
> 
> rm *tmp
> 
> # or add manually the column names: #SampleID	Run	Cruize	sample	marker	type 
> 
> ```
>
> 

Well done, now you can filter and make some ecological analysis within qiime or other plataform.

## Filtering count-data

save the follow script as `filter-features.sh`

how to use: `sbatch filter-features.sh table.gza 2 sampleMfile`

```bash
#!/bin/bash
#SBATCH --job-name=qiime
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=8
#SBATCH -t 6-00:00:00

# Load qiime

module load miniconda3-python-3.5
source activate quiime2-2018.0

# Remove low-frequency sequence

table=$1 # count-table processed with the denoise-paired part  
freq=$2
metadt=$3

prefix=$(ls $table | cut -d"_" -f1)

# 1 count-table filtering
qiime feature-table filter-features \
  --i-table $table \
  --p-min-frequency $freq \
  --o-filtered-table ${prefix}_table_${freq}_filtered.qza

# 2 rep-seqs filtering
qiime feature-table filter-seqs \
  --i-data ${prefix}_rep-seqs.qza \
  --i-table ${prefix}_table_${freq}_filtered.qza \
  --p-no-exclude-ids \
  --o-filtered-data ${prefix}_rep-seqs_table_${freq}_filtered.qza

exit
```

## Classify features

0) save this code with the name script `Classify-feature.sh`

> Example: sbatch Classify-feature.sh test_rep-seqs_table_3_filtered.qza 0.8 w2pr2_worms_API02.fasta.qza w2pr2_worms_API02.tax.qza

```bash
#!/bin/bash
#SBATCH --job-name=qiime
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# Load qiime

module load miniconda3-python-3.5
source activate quiime2-2018.0

query=$1 # feature-sequence 
identity=$2 # [from 0 to 1]
refreads=$3 # Unapped qza sequence reference
reftax=$4 # Taxonomy from sequence reference


tag=$(echo $reftax | cut -d'.' -f1)
taxonomy=$(echo $tag"_vs_"${query%.qza}.qza)

qiime feature-classifier classify-consensus-blast \
	--p-perc-identity $identity \
  --i-query $query \
  --i-reference-reads $refreads \
  --i-reference-taxonomy $reftax \
  --o-classification $taxonomy
  
exit
```

1) Remove min-frequency-features and exclude _unassign_, then visualize barplots

```bash
qiime feature-table summarize \
  --i-table  \
  --o-visualization w2pr2_worms_API02_vs_X04_rep-seqs_table_1_filtered.assigned.qza \
  --m-sample-metadata-file sampleMfile
  
table=X04_table_1_filtered.qza
seqs=X04_rep-seqs_table_1_filtered.qza
taxonomy=w2pr2_worms_API02_vs_X04_rep-seqs_table_1_filtered.qza  

table_filtered=${table%.qza}.assigned.qza  
seqs_filtered=${seqs%.qza}.assigned.qza 
tag=$(echo $seqs | cut -d'_' -f1)

qiime taxa filter-table \
  --i-table $table \
  --i-taxonomy $taxonomy \
  --p-exclude Unassigned \
  --o-filtered-table $table_filtered

qiime taxa filter-seqs \
  --i-sequences $seqs \
  --i-taxonomy $taxonomy \
  --p-exclude Unassigned \
  --o-filtered-sequences $seqs_filtered

# visualize

qiime taxa barplot \
  --i-table $table_filtered \
  --i-taxonomy $taxonomy \
  --m-metadata-file sampleMfile \
  --o-visualization ${tag}_filtered_assigned_barplot.qzv

exit
```

## Exporting Files

#### Convert taxonomy

>  Take less than 1-2 minute

We import these data into QIIME 2 Artifacts. Ex. the  reference taxonomy file (`w2pr2_worms_API02.tax`) is a tab-separated (TSV) file without a header, we must specify `HeaderlessTSVTaxonomyFormat` as the *source format* since the default *source format* requires a header.

1) In case of gapped fasta format, make a non gapped file. Example:

```bash
file=w2pr2_worms_API02.fasta

sed 'N;s/\n/\t/' $file | awk '{gsub(/[-, .]/, "", $2); print}' | awk '{print $1"\n"$2}' > w2pr2_worms_API02_no_gaps.fasta

# sanity check

grep '^>' w2pr2_worms_API02.fasta | sed 's/^>//' > ids.log
grep '^>' w2pr2_worms_API02_no_gaps.fasta | sed 's/^>//' > no_gaps_ids.log

diff ids.log no_gaps_ids.log

```

2) 

```bash
qiime tools import \
  --input-path w2pr2_worms_API02_no_gaps.fasta \
  --output-path w2pr2_worms_API02_no_gaps.fasta.qza \
  --type 'FeatureData[Sequence]'
  
# FeatureData[Sequence]
  
# 2 TSVTaxonomyFormat file

qiime tools import --type 'FeatureData[Taxonomy]' --input-format 'HeaderlessTSVTaxonomyFormat' --input-path w2pr2_worms_API02.tax --output-path w2pr2_worms_API02.tax.qza
	
```

Additional formats and import types:

```bash
qiime tools import --show-importable-types
qiime tools import --show-importable-formats
```

#### Convert outputs

> Example: `sbatch export-features.sh test_table_1_filtered.qza test_rep-seqs_table_1_filtered.qza w2pr2_worms_API02_vs_test_rep-seqs_table_1_filtered.qza` 

```bash
#!/bin/bash
#SBATCH --job-name=qiime
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=4
#SBATCH -t 6-00:00:00

# Load qiime

module load miniconda3-python-3.5
source activate quiime2-2018.0

count_feature=$1
seqs_feature=$2
classify_feature=$3


mkdir -p tabular_files

# Export feature-table
qiime tools export \
  --input-path $count_feature \
  --output-path tabular_files

biom convert \
	-i tabular_files/feature-table.biom \
	-o tabular_files/${count_feature%.qza}.tsv --to-tsv

# Export sequence table
qiime tools export \
  --input-path $seqs_feature \
  --output-path tabular_files

# Export Classification
qiime tools export \
	--input-path $classify_feature \
	--output-path tabular_files


exit

# ==== Visualizations
mkdir -p dataVis

# 1)
qiime feature-table summarize \
  --i-table ${prefix}_table.qza \
  --o-visualization dataVis/${prefix}_table.qzv \
  --m-sample-metadata-file $metadt

# 2)
qiime feature-table tabulate-seqs \
  --i-data ${prefix}_rep-seqs.qza \
  --o-visualization ${prefix}_rep-seqs.qzv

# 3)
qiime metadata tabulate \
  --m-input-file ${prefix}_denoising-stats.qza \
  --o-visualization ${prefix}_denoising-stats.qzv
  
  exit

# === within denoising script
# skipt this step after taxonomy part
mkdir -p tabular_files

# 3 export feature-table (biom)
qiime tools export \
  --input-path ${prefix}_table_${freq}_filtered.qza \
  --output-path tabular_files

# export sequence table (fasta)
qiime tools export \
  --input-path ${prefix}_rep-seqs_table_${freq}_filtered.qza \
  --output-path tabular_files

biom convert -i tabular_files/feature-table.biom -o tabular_files/${prefix}_feature_table_${freq}_filtered.tsv --to-tsv

# 4 visualize
qiime feature-table tabulate-seqs \
  --i-data ${prefix}_rep-seqs_table_${freq}_filtered.qza \
  --o-visualization ${prefix}_rep-seqs_table_${freq}_filtered.qzv

# 3.1 visualize
qiime feature-table summarize \
  --i-table ${prefix}_table_${freq}_filtered.qza \
  --o-visualization ${prefix}_table_${freq}_filtered.qzv \
  --m-sample-metadata-file $metadt

# within Classification step
# Export the taxononomy to tsv
mkdir -p tabular_files

qiime tools export \
	--input-path $taxonomy \
	--output-path tabular_files
	
exit

# Follow step is pipped to R
# Remove ====

# Remove Unassigned features from features table
qiime taxa filter-table \
  --i-table $table \
  --i-taxonomy $taxonomy \
  --p-exclude Unassigned \
  --o-filtered-table ${table%.qza}_assigned.qza

# Remove Unassigned features from sequences 
qiime taxa filter-seqs \
  --i-sequences $query \
  --i-taxonomy $taxonomy \
  --p-exclude Unassigned \
  --o-filtered-sequences ${query%.qza}_assigned.qza
  

# Export ======
mkdir -p tabular_files

# 1) Export the count-feature (first biom)

qiime tools export \
  --input-path ${table%.qza}_assigned.qza \
  --output-path tabular_files

biom convert -i tabular_files/feature-table.biom -o tabular_files/${table%.qza}_assigned.tsv --to-tsv

# 2) Export the taxononomy
qiime tools export \
	--input-path $taxonomy \
	--output-path tabular_files

# 3) Export fasta
qiime tools export \
  --input-path ${query%.qza}_assigned.qza \
  --output-path tabular_files
  

exit

```

workdir Tes:

```bash
/LUSTRE/bioinformatica_data/genomica_funcional/MG_18S/data/run04_20170418_18S/test
```

## Qiime complete pipeline

## XIXIM - 04

```bash

# 2) after visualize the f and r trimming positions run:
sbatch denoise-paired.sh 150 146 X04_demux_seqs.qza sampleMfile
sbatch denoise-paired.sh 150 149 X05_demux_seqs.qza sampleMfile
sbatch denoise-paired.sh 150 100 X06_demux_seqs.qza sampleMfile

# 3) Filter low-frequency features
sbatch filter-features.sh X04_table.qza 1 sampleMfile
sbatch filter-features.sh X05_table.qza 1 sampleMfile
sbatch filter-features.sh X06_table.qza 1 sampleMfile

# 4) Classifying features
sbatch Classify-feature.sh X04_rep-seqs_table_1_filtered.qza 0.8 w2pr2_worms_API02_no_gaps.fasta.qza w2pr2_worms_API02.tax.qza
#
sbatch Classify-feature.sh X05_rep-seqs_table_1_filtered.qza 0.8 w2pr2_worms_API02_no_gaps.fasta.qza w2pr2_worms_API02.tax.qza
#
sbatch Classify-feature.sh X06_rep-seqs_table_1_filtered.qza 0.8 w2pr2_worms_API02_no_gaps.fasta.qza w2pr2_worms_API02.tax.qza





```