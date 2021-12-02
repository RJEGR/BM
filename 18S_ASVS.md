## Trimming (amplicon-) primers from both ends of paired-end reads

[cite](https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads)

```bash
#wd=
python3 -m pip install --user --upgrade cutadapt
cd $wd
ln -s ~/.local/bin/cutadapt .
```

Configure it and save in a script of name `clip.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
 
for file in $(ls *R1*gz | grep -e 'fastq' -e 'fq' -e 'gz')
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fastq.gz}"
./cutadapt -a ^TTGTACACACCGCCC...GTAGGTGAACCTGCRGAAGG -A ^CCTTCYGCAGGTTCACCTAC...GGGCGGTGTGTACAA --discard-untrimmed --cores 24 -o ${base}_R1.clipped.fastq.gz -p ${base}_R2.clipped.fastq.gz ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz 
done

mkdir -p clipped
mv *clipped.fastq.gz ./clipped
exit 0
#cutadapt -a TTGTACACACCGCCC...GTAGGTGAACCTGCRGAAGG -A CCTTCYGCAGGTTCACCTAC...GGGCGGTGTGTACAA --discard-untrimmed -o Nombre_archivo_F -p Nombre_archivo_R input_file_F input_file_R

```

Check results

```bash
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

mkdir -p fastqc
fastqc *.gz -t 24 --nogroup -o ./fastqc 

exit 0
#export PATH=/LUSTRE/apps/Anaconda/conda2/bin:$PATH
#source activate multiqc_py2.7
#multiqc ./fastqc/*zip -o ./multiqc
```

count the number of reads

```bash
for file in $(ls *R1*gz | grep -e 'fastq' -e 'fq' -e 'gz')
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fastq.gz}"
#zcat ${base}_R1_001.fastq.gz | grep -c "^@M" ;
#zcat ${base}_R1.clipped.fastq.gz | grep -c "^@M" ; 
echo $base | awk 'BEGIN { FS = "-" } ; { print $3"."$2 }'
done

# 

for file in $(ls *R1*gz | grep -e 'fastq' -e 'fq' -e 'gz')
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fastq.gz}"
zcat ${base}_R1_001.qtrim.gz | grep -c "^@M" ;
#zcat ${base}_R1.clipped.fastq.gz | grep -c "^@M" ; 
echo $base | awk 'BEGIN { FS = "-" } ; { print $3"."$2 }'
done


#

for file in $(ls *R1*gz | grep -e 'fastq' -e 'fq' -e 'gz'); do withpath="${file}"; filename=${withpath##*/}; base="${filename%*_R*.qtrim.gz}"; echo $base; done
```



slurm the script `dada.slurm dada workdir 100 80 `

```bash
#!/bin/bash
#SBATCH --job-name=dada_18S
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"


TESTNAME=$1
LIBRARIES_DIR=$2
trunF=$3
trunR=$4

module load R-3.5.0
Rscript --vanilla multirun_18S.R $TESTNAME $LIBRARIES_DIR $trunF $trunR $SLURM_NPROCS 

exit 0


```

