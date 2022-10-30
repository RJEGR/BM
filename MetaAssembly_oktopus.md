Add: Lobulo optico de hembras y Corazon Macho.



Work directory in:

```bash
/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oktopus_full_assembly
```



- 118 archivos
  - 59 librerias (pareadas)



0) Run the filtering

```bash
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=qtrim 
#SBATCH --output=trm-%j.log 
#SBATCH --error=trm-%j.err 
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24 
#SBATCH -t 06-00:00:00 
#SBATCH --exclusive

TRIMMOMATIC=/LUSTRE/bioinformatica_data/RNA/ricardo/bioinformatics/Trimmomatic-0.36
TRUSEQ=/home/rgomez

for file in $(ls *R1*gz | grep -e 'fastq' -e 'fq') 
do 
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fq.gz}" 
java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 -threads $SLURM_NPROCS \
    ${base}_R1.fq.gz ${base}_R2.fq.gz \
    ${base}_R1.P.qtrim.fq.gz ${base}_R1.UP.qtrim.fq.gz \
    ${base}_R2.P.qtrim.fq.gz ${base}_R2.UP.qtrim.fq.gz \
    ILLUMINACLIP:$TRUSEQ/TruSeq3-PE-2.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:36 LEADING:5 TRAILING:5
done
```

1) Prepare meta-data

```bash
ls *R2* | sort > R2.tmp && ls *R1* | sort > R1.tmp && \
cut -d "-" -f1 R1.tmp > factors.tmp && \
cut -d "_" -f1 R1.tmp > codes.tmp && \
paste factors.tmp codes.tmp R1.tmp R2.tmp | awk '{gsub(/\-/,"_",$2); \
print $1,$2,$3,$4}' | column -t > samples.file && rm *tmp

# Check data
wc -l samples.file

```

2) Run the assembly:

**Hay un problema con los files dentro de ASSEMBLY_COMPLETE**

for i in $(ls *gz); do unlink $i; done



ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/OvG/*.gz .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/OG_Female/*.gz .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/OG_Female/*.gz .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/WB_male/*.gz .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/OL_Male/*gz .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/TESTIS/*gz .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/WB_female/*gz .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/PULPO_Project/Heart_Male_relabeled/* .

ln -s /LUSTRE/bioinformatica_data/genomica_funcional/diana/secuencias/* .



rm slurm*

sbatch ./assembly.sh -f samples.file -m 100 -p 24 -o full



```bash
--grid_exec <string>                 :your command-line utility for submitting jobs to the grid.
# This should be a command line tool that accepts a single parameter:
# ${your_submission_tool} /path/to/file/containing/commands.txt
# and this submission tool should exit(0) upon successful 
# completion of all commands.
#
--grid_node_CPU <int>                number of threads for each parallel process to leverage. (default: 1)
#
--grid_node_max_memory <string>         max memory targeted for each grid node. (default: 1G)
#
#            The --grid_node_CPU and --grid_node_max_memory are applied as 
#              the --CPU and --max_memory parameters for the Trinity jobs run in 
#              Trinity Phase 2 (assembly of read clusters)
#
#

export PATH=/LUSTRE/bioinformatica_data/biocomp/gmelendrez/tools/bin/HpcGridRunner:$PATH

grid_conf=/LUSTRE/bioinformatica_data/biocomp/gmelendrez/tools/bin/HpcGridRunner/hpc_conf/SLURM.test.conf

grid_path=/LUSTRE/bioinformatica_data/biocomp/gmelendrez/tools/bin/HpcGridRunner/

Trinity ...  --grid_exec "$grid_path/hpc_cmds_GridRunner.pl --grid_conf $grid_conf -c"



```



```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=assembly
#SBATCH --output=slrm-%j.log
#SBATCH --error=slrm-%j.err
#SBATCH -N 3
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -p cicese

# ======================================================
# Step 1: in silico normalization step and kmer counting
# gcc-7.2 needed for insilico_read_normalization

wd=$1

module load gcc-7.2

export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/:$PATH
export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/trinity-plugins/jellyfish-2.2.6/bin/:$PATH

which Trinity

Trinity --seqType fq --max_memory 100G --samples_file samples.file --no_normalize_reads --CPU 24 --output trinity_out --no_salmon

# Trinity --seqType fq \
         --left R1.fastq  \
         --right R2.fastq\
         --max_memory 100G \
         --output trinity_${wd} \
         --CPU $SLURM_NPROCS \
         --no_run_inchworm


# ===================
# Step 2: run inchworm, stop before chrysalis:

wd=$(echo trinity_*)


export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/Inchworm/bin/:$PATH
export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/support_scripts/:$PATH
export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/Chrysalis/bin/:$PATH

which inchworm

inchworm --kmers $wd/jellyfish.kmers.fa \
				 --run_inchworm -K 25 -L 25 \
				 --monitor 1   --DS  \
				 --num_threads $SLURM_NPROCS \
				 --PARALLEL_IWORM  > $wd/inchworm.K25.L25.DS.fa.tmp

mv $wd/inchworm.K25.L25.DS.fa.tmp $wd/inchworm.K25.L25.DS.fa

touch $wd/inchworm.K25.L25.DS.fa.finished


# support_scripts
# filter_iworm_by_min_length_or_cov.pl IwormFastaFile min_length min_cov
which filter_iworm_by_min_length_or_cov.pl

filter_iworm_by_min_length_or_cov.pl $wd/inchworm.K25.L25.DS.fa 100 10 > $wd/chrysalis/inchworm.K25.L25.DS.fa.min100

bowtie2-build -o 3 --threads $SLURM_NPROCS $wd/chrysalis/inchworm.K25.L25.DS.fa.min100 $wd/chrysalis/inchworm.K25.L25.DS.fa.min100 1>/dev/null

bash -c " set -o pipefail;bowtie2 --local -k 2 --threads $SLURM_NPROCS -f --score-min G,46,0 -x $wd/chrysalis/inchworm.K25.L25.DS.fa.min100 $wd/both.fa | samtools view -@ $SLURM_NPROCS -F4 -Sb - | samtools sort -m 2236962133 -@ $SLURM_NPROCS -no - - > $wd/chrysalis/iworm.bowtie.nameSorted.bam"


# support_scripts
which scaffold_iworm_contigs.pl

scaffold_iworm_contigs.pl $wd/chrysalis/iworm.bowtie.nameSorted.bam $wd/inchworm.K25.L25.DS.fa > $wd/chrysalis/iworm_scaffolds.txt

# Chrysalis
# GraphFromFasta: Makes a graph out of a fasta

which GraphFromFasta

GraphFromFasta -i $wd/inchworm.K25.L25.DS.fa -r $wd/both.fa -min_contig_length 200 -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05 -t 24 -k 24 -kk 48  -scaffolding $wd/chrysalis/iworm_scaffolds.txt  > $wd/chrysalis/iworm_cluster_welds_graph.txt

which BubbleUpClustering

BubbleUpClustering -i $wd/inchworm.K25.L25.DS.fa  -weld_graph $wd/chrysalis/iworm_cluster_welds_graph.txt -min_contig_length 200  > $wd/chrysalis/GraphFromIwormFasta.out

exit

```



## Exploratory data analysis 

- <u>Step 0</u>: **in silico normalization step and kmer counting:**  (Trinity Parameters:   `--seqType fq --left R1.fastq --right R2.fastq --max_memory 100G --CPU 48 —no_run_inchworm`)
  - `$wd/jellyfish.kmers.fa.histo`
  - `$wd/both.fa.read_count`

*Inchworm* assembles the RNA-seq data into the unique sequences of transcripts, often generating full-length transcripts for a dominant isoform, but then reports just the unique portions of alternatively spliced transcripts.

- <u>Step 1</u>: **run inchworm**, stop before chrysalis:
  - `$wd/chrysalis/iworm.bowtie.nameSorted.bam`
  - `$wd/chrysalis/iworm_cluster_welds_graph.txt`

*Chrysalis* clusters the Inchworm contigs into clusters and constructs complete de Bruijn graphs for each cluster. Each cluster represents the full transcriptonal complexity for a given gene (or sets of genes that share sequences in common). Chrysalis then partitions the full read set among these disjoint graphs.

- Step 2: Run Chrystalis



*Butterfly* then processes the individual graphs in parallel, tracing the paths that reads and pairs of reads take within the graph, ultimately reporting full-length transcripts for alternatively spliced isoforms, and teasing apart transcripts that corresponds to paralogous genes.

- Step 3: 
  - kj



