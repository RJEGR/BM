Ingresar con tus credenciales:

Ejemplo:

`ssh usuario@omica`

Y vallamos a nuestro directorio de trabajo:

`cd /LUSTRE/bioinformatica_data/genomica_funcional/RAD2020`

## Para Stacks

Estamos muy emocionados por nuestros datos de genomica poblacional. Sin embargo nos encontramos con el reto de que nuestros analisis bioinformaticos demoran mas de lo esperado, para ello usaremos un servidor de computo (Clusters) que tiene la capacidad de hacer analisis complejos en un menor tiempo. 

`Slurm`, es un sistema de gestión de tareas que nos permite hacer trabajos dentro de un clusters de computo de manera optima. Este sistema de tareas viene contenido dentro del `omica` en el cual ingresamos previamente.

Creamos un script con las directrireces necesarias ejecutar nuestros comandos atraves del administrador de tareas. Este script lo guardamos en un archivo de texto plano llamado `process_radtags.sh` y ejecutamos nuestra tarea (`sbatch`) como a continuación:

`sbatch process_radtags.sh raw/UO_C716_1.fastq.gz barcodes/barcode_C716_AM.txt my_out_folder`

El script `process_radtags.sh` va a solicitar tres variables de entrada:

1. Archivo fasta (Ej.`raw/UO_C716_1.fastq.gz`)
2. Lista de barcodes (Ej. `barcodes/barcode_C716_AM.txt`)
3. nombre de archivo de salida (Ej. `my_out_folder`)

No sera necesario sobre-escribir dentro del script las variables mencioadas debido a que estas seran remplazadas por lo que se describe dentro de la sintaxis del sbatch.

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



### Otras notas

Como descomprimimos el archivo de trabajo

```bash
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=tar 
#SBATCH --output=trm-%j.log 
#SBATCH --error=trm-%j.err 
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24 
#SBATCH --exclusive

tar -zxvf FLafarga_CICESE_Abalone_20161108−01479.tar.gz raw

exit

# tar -zxvf config.tar.gz etc/default/sysstat
```



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



