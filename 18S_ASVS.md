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
#SBATCH -N 3
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

for file in $(ls *R1*gz | grep -e 'fastq' -e 'fq' -e 'gz')
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fastq.gz}"
./cutadapt -a ^TTGTACACACCGCCC...GTAGGTGAACCTGCRGAAGG -A ^CCTTCYGCAGGTTCACCTAC...GGGCGGTGTGTACAA --discard-untrimmed -o ${base}_R1.clipped.fastq.gz -p ${base}_R2.clipped.fastq.gz ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz
done

mkdir -p clipped
mv *clipped.fastq.gz ./clipped
exit 0
#cutadapt -a TTGTACACACCGCCC...GTAGGTGAACCTGCRGAAGG -A CCTTCYGCAGGTTCACCTAC...GGGCGGTGTGTACAA --discard-untrimmed -o Nombre_archivo_F -p Nombre_archivo_R input_file_F input_file_R

```

#