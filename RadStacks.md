Ingresar con credenciales

En esta ocasión el proceso de instalación no presentó ninguna complicación, ya quedó la versión 2.5 de stacks, se encuentra en la trayectoria:

```bash
# path
# which process_radtags
# ls /LUSTRE/apps/bioinformatica/stacks-2.5/bin 

export PATH=/LUSTRE/apps/bioinformatica/stacks-2.5/bin/:$PATH

module load gcc-7.2
```



| N    | A.             |
| ---- | -------------- |
| 14   | Amarillo       |
| 16   | Rojo (Cultivo) |
| 20   | Azul           |
| 16   | Negro          |
|      |                |
|      |                |

```bash
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=fqc 
#SBATCH --output=trm-%j.log 
#SBATCH --error=trm-%j.err 
#SBATCH -N 3
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24 
#SBATCH -t 06-00:00:00

fasta=$1
outdir=${fasta%.fastq.gz}_fastqc

mkdir $outdir

fastqc $fasta -t $SLURM_NPROCS --nogroup -o $outdir
```

`sbatch fastqc.sh UO_C716_1.fastq.gz`

```bash
mkdir multiqc
export PATH=/LUSTRE/apps/Anaconda/conda2/bin:$PATH
source activate multiqc_py2.7
multiqc ./fastqc/*zip -o ./multiqc --data-format json --export
```



for stacks

```bash
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=rtgs 
#SBATCH --output=trm-%j.log 
#SBATCH --error=trm-%j.err 
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24 
#SBATCH -t 06-00:00:00 
#SBATCH --exclusive

which process_radtags

export PATH=/LUSTRE/apps/bioinformatica/stacks-2.5/bin/:$PATH

module load gcc-7.2

fasta=$1 # Ex. UO_C716_1.fastq.gz
bars=$2 # Ex. barcode_C716_AM.txt

outdir=$3 # Ex. my_out_folder

mkdir $outdir

process_radtags -f $fasta -b $bars -o $outdir -e 'sbfI' -c -q -r -t 75

exit

```

`sbatch process_radtags.sh raw/UO_C716_1.fastq.gz barcodes/barcode_C716_AM.txt my_out_folder`



