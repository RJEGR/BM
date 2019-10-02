# Análisis de datos de secuenciación masiva de pequeños RNAs no-codificantes

>  27 - 28 de septiembre en CICESE
>
> PhD Cei Abreu, LANGEBIO

> Descripción:

Todos los seres vivos, además de transcribir genes que codifican para proteínas, lo hacen para
una amplia variedad de genes llamados no-codificantes. Muchos organismos también producen
RNAs a partir de templados de RNA. Por razones prácticas los RNAs no-codificantes se
clasifican en cortos (pequeños) y largos. Quizá la clase más famosa de RNAs pequeños son
los microRNAs, por tratarse de moléculas esenciales para el desarrollo de plantas y animales.
La tecnología de secuenciación masiva más popular de los últimos años (Illumina) permite
exploraciones muy completas de los RNAs pequeños de una amplia diversidad de muestras,
desde organismos u órganos completos, material extracelular, organismos en interacción o
inclusive células individuales. En este curso los participantes aprenderán a apreciar las
peculiaridades y limitaciones de la secuenciación de RNAs pequeños. Usando R dentro de
RStudio, y algunas herramientas en Linux, aprenderán también a analizar y explorar datos de
secuenciación masiva de RNAs pequeños, así como interpretar los resultados obtenidos.

## Instrucciones

#### Ingresando al cluster

0) Desde windows ingresa a la herramienta **PuTTY**, en mac y linux desde tu terminal

1) Selecciona una cuenta e ingresa al cluster

```bash
ssh curso01@omica
password: RNA-cicese01
```

| Cuenta  |    Clave     |
| :-----: | :----------: |
| curso01 | RNA-cicese01 |
| curso02 | RNA-cicese02 |
| curso03 | RNA-cicese03 |
| curso04 | RNA-cicese04 |
| curso05 | RNA-cicese05 |
| curso06 | RNA-cicese06 |
| curso07 | RNA-cicese07 |
| curso08 | RNA-cicese08 |
| curso09 | RNA-cicese09 |
| curso10 | RNA-cicese10 |
| curso11 | RNA-cicese11 |
| curso12 | RNA-cicese12 |
| curso13 | RNA-cicese13 |
| curso14 | RNA-cicese14 |
| curso15 | RNA-cicese15 |
| curso16 | RNA-cicese16 |
| curso17 | RNA-cicese17 |
| curso18 | RNA-cicese18 |
| curso19 | RNA-cicese19 |
| curso20 | RNA-cicese20 |
| curso21 | RNA-cicese21 |
| curso22 | RNA-cicese22 |
| curso23 | RNA-cicese23 |
| curso24 | RNA-cicese24 |
|         |              |



2) Una vez ingresado al cluster, ingresamos las siguientes lineas de comando:

```
bash ~/accesa-nodo
cd curso2019
```

3) Revisamos los archivos con los que trabajaremos:

```bash
ls
```

4) Realizamos la copia del directorio `data_min`, **en caso de no encontrar la carpeta dentro de tu directorio de trabajo**:

```bash
cp -r /LUSTRE/bioinformatica_data/curso2019/accounts/cei/data_min .
cd data_min
```

5) Revisamos que tengamos la ruta de progamas en nuestra sesión `echo $PATH` . **En caso de no contar con la ruta de los programas con los que trabajaremos:**

```bash
export PATH=/LUSTRE/bioinformatica_data/genomica_funcional/cei/bin$PATH
```

6) Vamos al tutorial del curso en:

[Prácticas de Bioinformática para el análisis de datos de sRNA-Seq](http://datos.langebio.cinvestav.mx/~cei/cursos/CICESE/)



# Trabajando con R (version 3.5 o superior)

**0) Instalación** 

```R
install.packages("BiocManager")

BiocManager::install()

.pkgs <- c("Rsamtools", "GenomicFeatures", "GenomicAlignments",
        "ggplot2", "pheatmap", "RColorBrewer","AnnotationDbi",
        "rtracklayer","GenomicRanges","edgeR","statmod","readr","biomaRt",
        "Biostrings","Biobase","knitr","affy","RCurl")

BiocManager::install(.pkgs)
```

**1) Cargando los paquetes**

```R
.inst <- .pkgs %in% installed.packages()

if(any(!.inst)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(.pkgs[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.pkgs), require, character.only = TRUE)
```





## **Linux tools** 

- samtools [http://www.htslib.org/](http://www.htslib.org/)
- bowtie1 [http://bowtie-bio.sourceforge.net/index.shtml](http://bowtie-bio.sourceforge.net/index.shtml)
- ShortStack https://github.com/MikeAxtell/ShortStack



```bash
# Para usarlo agregar esta trayectoria a la variable ambiental PATH:
export PATH=/LUSTRE/apps/bioinformatica/Shortstack:$PATH
Shortstack

#Este requiere de las otras aplicaciones que solicitaban, bowtie1 y samtools, del primero se encuentran las versiones 1.1.2 y 1.2.1.1 en las carpetas:
# <-- versión 1.1.2
export PATH=/LUSTRE/apps/bioinformatica/bowtie1:$PATH   

# <-- versión 1.2.1.1
export PATH=/LUSTRE/apps/bioinformatica/bowtie-1.2.1.1:$PATH  

# y samtools
export PATH=/LUSTRE/apps/bioinformatica/samtools-1.7/bin:$PATH

# Para usarlos en conjunto con Shortstack

export PATH=/LUSTRE/apps/bioinformatica/Shortstack:/LUSTRE/apps/bioinformatica/samtools-1.7/bin:/LUSTRE/apps/bioinformatica/bowtie-1.2.1.1:$PATH`
```

