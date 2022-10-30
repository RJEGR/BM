### Metaensamble de pulpo

Ricardo Gómez Reyes

**Reto:**

Tuve problemas al ejecutar transrate dentro del cluster para analizar el transcriptoma completo de pulpo. En principio esto se debe a la cantidad de datos que estamos usando para elaborar el meta-ensamble (estamos hablando de más de 200GB de datos que se ensamblan), apensar de modificar los parametros de transrate para sugerir una mayor cantidad de memoria, la memoria colapsa y se aborta el análisis. Así que, no se está logrando filtrar los contigs del ensamble con dicha herramienta. Por consiguiente se filtraron contigs basado en el criterio de la cobertura (Para evaluar la composición de nuestro ensamble, queremos capturar y contar todas las lecturas paired-end que se asignan a los andamios ensamblados, incluidas las correctamente "orientadas" y las que no), ejecutamos el proceso a continuación:  

```bash
# Previamente concatenamos nuestras bibliotecas F y R en un solo archivo *fq

# 1
bowtie2-build Trinity.fasta Trinity.fasta

# 2
srun bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -1 LOP_H_F.fq -2 LOP_H_R.fq | samtools view -@10 -Sb -o bowtie2.bam

# 3

bamtools stats -in bowtie2.bam

```

- Los resultados de la prueba fueron los siguientes:

Total reads:    69,365,541
**Mapped reads:    69,365,541 (100%)**
**Both pairs mapped: 54,087,779 (77.975%)**

Forward strand:   33,138,976 (47.7744%)
Reverse strand:   36,226,565 (52.2256%)

Adicionalmente, se estimo la abundancia de los transcritos evaluados en `properly_mapped_pairs` y se filtro en base a la abundancia (tpm) > 1. Los resultados en comparacion con el ensamble crudo se presentan a continuación

|                            | Trinity.fasta        | properly_mapped_pairs                         |
| -------------------------- | -------------------- | --------------------------------------------- |
| Total trinity 'genes':     | 384,151.00           | 99,853 (26%, con respecto al ensamble crudo)  |
| Total trinity transcripts: | 895,197.00           | 361,343 (40%, con respecto al ensamble crudo) |
| Percent GC:                | 37.49                | 37.09                                         |
| Contig N10:                | 2266                 | 2860                                          |
| Contig N20:                | 1496                 | 2016                                          |
| Contig N30:                | 1067                 | 1553                                          |
| Contig N40:                | 785                  | 1219                                          |
| Abundancia > 1:            | 78.32 % (701,083.00) | 85.69% (309,633.00)                           |
|                            |                      |                                               |

Transrate usa ambas métricas para evaluar calidad de los contigs cobertura y abundancia. Si se consiera que el porcentaje de contigs que se suelen filtrar con *transrate* pueden ir del 40 % al 60% de todo el ensamble, los resultados obtenidos por el método en cuestión son de esperarse.

Con esta prueba se puede considerar usar los contigs que fueron mapeados por la información pareada de lecturas en su orientación correcta (*Properly_mapped_pairs*), de esta manera mitigamos quimerismo, inserciones y redundancia dentro de nuestro ensamble. 

|        Clean abundance > 1 |               |                       |
| -------------------------: | ------------- | --------------------- |
|                            | Trinity.fasta | properly_mapped_pairs |
|     Total trinity 'genes': | 339262        | 96849                 |
| Total trinity transcripts: | 701083        | 309633                |
|        % of loose assembly | 21.68%        | 14.31%                |
|                Percent GC: | 37.6          | 37.04                 |
|                Contig N10: | 2083          | 2654                  |
|                Contig N20: | 1358          | 1870                  |
|                Contig N30: | 946           | 1434                  |
|                Contig N40: | 693           | 1119                  |
|                Contig N50: | 528           | 877                   |
|                            |               |                       |
|       Median contig length | 309           | 463                   |
|            Average contig: | 465.52        | 675.24                |
|     Total assembled bases: | 326,366,757   | 209,076,316           |



El problema con transrate

```
Se modificó la configuración de snap (snap.rb); agregando un numero mayor en el parametro mcp (cmd << " -mcp 10000000"); El archivo modificado se encuentra en: /LUSTRE/bioinformatica_data/RNA/ricardo/bioinformatics/transrate-1.0.3-linux-x86_64/lib/app/lib/transrate/snap.rb sin éxito en la ejecución de transrate con ensamble de lecturas abundantes.


```

