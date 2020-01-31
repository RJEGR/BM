# Stacks: Construyendo loci a partir de secuencias de lectura corta

> 29-30 de Enero, 2020 en CICESE
>
> Ponente: Dr. Cristián Araneda, Laboratorio de Genética y Biotecnología en Acuicultura, Universidad de Chile

## Para empezar

Ingresar con tus credenciales:

`ssh usuario@omica`

`Password:******`

Y vallamos a nuestro directorio de trabajo:

`cd /LUSTRE/bioinformatica_data/genomica_funcional/RAD2020`

Vamos a cargar los programas que usaremos. Versión 2.5 de stacks que se encuentra en la trayectoria:

```bash
export PATH=/LUSTRE/apps/bioinformatica/stacks-2.5/bin/:$PATH
module load gcc-7.2
```

Para dejar de manera permanente el programa alojado en nuestras variables de ambiente usemos el siguiente comando. **De esta manera no sera necesario exportar la ruta del programa cada vez que iniciemos sesion dentro del cluster:**

`echo 'export PATH=/LUSTRE/apps/bioinformatica/stacks-2.5/bin/:$PATH' >> ~/.bash_profile`

**Aviso de privacidad: Ningún dato presente en este documento compromete resultados de investigación privados**

## 1. Procesando Radtags

Estamos muy emocionados por nuestros datos de genomica poblacional. Sin embargo nos encontramos con el reto de que nuestros analisis bioinformaticos demoran mas de lo esperado, para ello usaremos un servidor de computo (Clusters) que tiene la capacidad de hacer analisis complejos en un menor tiempo. 

`Slurm`, es un sistema de gestión de tareas que nos permite hacer trabajos dentro de un clusters de computo de manera optima. Este sistema de tareas viene contenido dentro del `omica` en el cual ingresamos previamente.

Creamos un script con las directrireces necesarias ejecutar nuestros comandos atraves del administrador de tareas. Este script lo guardamos en un archivo de texto plano llamado `process_radtags.sh` y ejecutamos nuestra tarea (`sbatch`) como a continuación:

`sbatch process_radtags.sh raw/UO_C716_1.fastq.gz barcodes/barcode_C716_AM.txt my_out_folder 100` 

El script `process_radtags.sh` va a solicitar tres variables de entrada:

1. Archivo fasta (Ej.`raw/UO_C716_1.fastq.gz`)
2. Lista de barcodes (Ej. `barcodes/barcode_C716_AM.txt`)
3. nombre de archivo de salida (Ej. `my_out_folder`)
4. Truncate final read length to this value (Ej. `75-100`)

No sera necesario sobre-escribir dentro del script las variables mencioadas debido a que estas seran remplazadas por lo que se describe dentro de la sintaxis del sbatch.

```bash
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=rtgs 
#SBATCH --output=rtgs-%j.log 
#SBATCH --error=rtgs-%j.err 
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24 
#SBATCH -t 06-00:00:00 
#SBATCH --exclusive

export PATH=/LUSTRE/apps/bioinformatica/stacks-2.5/bin/:$PATH

module load gcc-7.2
which process_radtags

fasta=$1 # Ex. UO_C716_1.fastq.gz
bars=$2 # Ex. barcode_C716_AM.txt

outdir=$3 # Ex. my_out_folder

len=$4

mkdir -p $outdir

process_radtags -f $fasta -b $bars -o $outdir -e 'sbfI' -c -q -r -t $len

exit
```



## 2. Ustacks (Construccion de SNPs)

1. **Algunas pruebas previas**: vamos a variar el minimo de covertura (m), maximo de distancia entre stacks (M) y  maximo de distancia (N)  para evaluar el procesamiento de stacks que se generen hasta que encontremos nuestros parametros de _equilibrio_ que mas se acerque a la columna de distribución normal entre los valores de las pruebas. Vamos a usar aquella muestra cuyo numero de lecturas representen al promedio de las lecturas de todas las muestras procesadas con `process_radtags`. 1) Vamos a evaluar el parametro M hasta tomar el valor normal 2) Variamos los valores n/M.  Aqui tenemos que ser cuidadonos porque nuestro numero de stacks esperados tiene que ser consistente con las medidas realizadas previamente con los calculos de contenido GC genomico y corte enzimatico (Por ejemplo usando la herramienta `simRAD` en R).

| 8    | ..   | 1    | M    | # Staks |
| ---- | ---- | ---- | ---- | ------- |
| 10   | ..   | 5    | N    | # Staks |
| 7    | ..   | 2    | m    | # Staks |

Creamos un script con las directrireces necesarias ejecutar nuestros comandos atraves del administrador de tareas. Este script lo guardamos en un archivo de texto plano llamado `ustacks.sh` y ejecutamos nuestra tarea (`sbatch`) como a continuación:

`sbatch ustacks.sh /path/to/file.fq.gz 9 7 2 4` 

```bash
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=ustks
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

# How to:
# sbatch ustacks.sh /path/to/file.fq.gz 9 7 2 4

export PATH=/LUSTRE/apps/bioinformatica/stacks-2.5/bin/:$PATH

module load gcc-7.2

which ustacks

# i : a unique integer ID for this sample.
# -m — Minimum depth of coverage required to create a stack (default 3).
# -M — Maximum distance (in nucleotides) allowed between stacks (default 2).
# -N — Maximum distance allowed to align secondary reads to primary stacks (default: M + 2).
# --alpha [num] — chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001 
# d --deleverage: enable the Deleveraging algorithm, used for resolving over merged tags.
# --model_type 'snp' (default)

fst=$1
i=$2
m=$3 # Ej. 7
M=$4 # Ej. 2
N=$5 # Ej. 4

mkdir -p $PWD/m${m}M${M}N${N}n${n}_snp

ustacks -t gzfastq -f $fst  -o $PWD/m${m}M${M}N${N}n${n}_snp -i $i -m $m -M $M -N $N --alpha 0.05 -d -H -p $SLURM_NPROCS

exit
```

Para esta prueba usaremos una muestra representativa: `/LUSTRE/bioinformatica_data/genomica_funcional/RAD2020/samples_an/AN_9.fq.gz`

Para nuestro datos, notamos que 6 fue el valor normal para m y 4,6 los valores esperados para M y N. Posteriormente corremos usando estos parametros estándares el modulo ustacks para el resto de nuestros datos:

2. Estandarizar y ejecutar recursivamente a todos las bibliotecas

Vamos a re-muestrear cada biblioteca *para asegurar la misma contribución de variacion en la secuencia de todas las muestras*. Esto nos dara la misma profundidad en numero de lecturas para todas nuestras muestras antes de procesar nuestros `stacks`

```bash
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=ustks
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24 
#SBATCH -t 06-00:00:00 

# Re-sampling
for i in *gz; do zcat $i | head -n10000000 > ${i%.fq.gz}.fq; done

ls *.fq > std_samples

```

Y finalmente ejecutamos `stacks` en todas las muestras:

```bash
#!/bin/bash

mkdir -p ustacks_gore

i=1

for file in $(cat std_samples)
do
    ustacks -p 24 -t fastq -m 6 -i $i -M 4 -N 6 \
            -f $file \
            -o ustacks_gore
    let "i+=1";
done
```



## 3. Construccion de catalogo

Ahora vamos a conserva el numero de alelos, SNPs y Loci dentro de un catalogo usando `cstacks`

```bash
cstacks -p 20 -o ustacks_gore -n 2 -s ustacks_gore/AN_4.fq.gz_2.5 -s ustacks_gore/AN_5.fq.gz_2.5 -s ustacks_gore/AN_6.fq.gz_2.5 -s ustacks_gore/AN_9.fq.gz_2.5
# Final catalog contains 56131 loci.

```

Finalmente genotipamos las muestras de un catalogo con `sstacks`

```bash
#!/bin/bash

mkdir -p ustacks_gore

i=1

for file in $(cat std_samples)
do
    sstacks -p 24 -c ustacks_gore \
            -s $file \
            -o ustacks_gore
    let "i+=1";
done
```

Para cada muestras vamos a encontrar archivos con extension `*matches.tsv` que corresponde a los SNPs que estan alineando con el catalogo.

## 4. Population



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

tar -xvf FLafarga_CICESE_Abalone_20161108−01479.tar.gz

exit

# tar -zxvf config.tar.gz etc/default/sysstat
```

Metricas de calidad:

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
multiqc ./*_fastqc/*zip -o ./multiqc --data-format json --export
```

> Notas elaboradas por Ricardo Gomez (Sin editar).

