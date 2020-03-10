## Coral DEGs

* fasta de DEGs a transdecoder: 
> Summary:
* Trinity.fasta: 226,049 sequence (`grep -c "^>"`)
* input : 766 DEGS.list ( `sort DEGS.list | uniq -c | wc -l `)
* 1449 longest_orfs.pep (`grep -c "^>" `)
  * of which 593 are unique 
* 788 predicted ORFs (`grep -c "^>" DEGS.Trinity.fasta.transdecoder.pep`)
  * of wichi 580 uniques 
* See experimental test to further analysis
```bash
# Rename headers with single name
awk '/^>/{print $1, $2; next}{print}' < Trinity.fasta | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > Trinity.renameHeaders.fasta

# 2 capturamos solo los DEGs

awk '{print $1}' DEGS.list | xargs -I {} grep -A1 "{}" Trinity.renameHeaders.fasta > DEGS.Trinity.fasta
```
```bash
 grep ">" DEGS.Trinity.fasta | awk '{print $1}' | sed 's/>//' > DEGS.fasta.list
 diff DEGS.fasta.list DEGS.list # 8 transcritos extras en el DEGS.Trinity.fasta ?? edit manually:

< TRINITY_DN39915_c0_g1_i10
< TRINITY_DN39915_c0_g1_i11
< TRINITY_DN39915_c0_g1_i12
< TRINITY_DN39915_c0_g1_i13
< TRINITY_DN39915_c0_g1_i14
< TRINITY_DN39915_c0_g1_i15
< TRINITY_DN40343_c2_g1_i10

```
Implementamos Transdecoder para predecir el marco de lectura:

```bash
DECODER=/LUSTRE/apps/bioinformatica/TransDecoder-3.0.1/
UTILS=/LUSTRE/apps/bioinformatica/trinityrnaseq/util/support_scripts
$UTILS/get_Trinity_gene_to_trans_map.pl DEGS.Trinity.fasta > DES.Trinity.gene_trans_map
#1 
srun $DECODER/TransDecoder.LongOrfs -G universal -t DEGS.Trinity.fasta --gene_trans_map DEGS.Trinity.gene_trans_map > DEGS_transdecoder.log &
#2

srun $DECODER/TransDecoder.Predict -t DEGS.Trinity.fasta --cpu 24 >> DEGS_transdecoder.log & 

# transdecoder is finished.  See output files DEGS.Trinity.fasta.transdecoder.*
```
## Experimental test (Including homology searches as ORF retention criteria): 

Para maximizar aún más la sensibilidad para capturar los ORF que pueden tener un significado funcional, se puede escanear todos los ORF en busca de homología con proteínas conocidas y retener todos esos ORF. Esto se puede hacer de dos formas populares: una búsqueda BLAST en una base de datos de proteínas conocidas y la búsqueda de PFAM para identificar dominios de proteínas comunes. En el contexto de TransDecoder, esto se hace de la siguiente manera:

Retener el mejor ORF por transcrito (`--single_best_orf`) basado en los resultados de dominios conservados evaluados con HMMER y la base de datos PFam (`retain_pfam_hits`) o blastp (`--retain_blastp_hits`) y la base de peptidos _uniprot_sprot_

```bash

cd trandecoder_dir

# HMMSCAN || Pfam-A.hmm
HMMSCAN=/LUSTRE/bioinformatica_data/genomica_funcional/bin/hmmer-3.1b2-linux-intel-x86_64/binaries
srun $HMMSCAN/hmmscan --cpu 24 --domtblout DEGsPFAM.out ../../Pfam-A.hmm longest_orfs.pep &

# BLASTP || uniprot_sprot
srun blastp -query longest_orfs.pep  \
    -db ../../uniprot_sprot.pep  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6 &
# retain_pfam_hits and single_best_orf:
srun 
$DECODER/TransDecoder.Predict -t DEGS.Trinity.fasta --retain_pfam_hits DEGS.Trinity.fasta.transdecoder_dir/DEGsPFAM.out --single_best_orf

mv DEGS.Trinity.fasta.transdecoder.pep DEGS.Trinity.fasta.transdecoder.pfamHits.pep

# or --retain_blastp_hits and and single_best_orf

$DECODER/TransDecoder.Predict -t DEGS.Trinity.fasta --retain_blastp_hits DEGS.Trinity.fasta.transdecoder_dir/blastp.outfmt6 --single_best_orf

mv DEGS.Trinity.fasta.transdecoder.pep DEGS.Trinity.fasta.transdecoder.
blastpHits.pep
```

**Resultados**

* 663 ORFs totales (`grep -c ">" DEGS.Trinity.fasta.transdecoder.pfamHits.pep`) usando PFAM
  * 582 unicos
* 662 ORFs  totales (`grep -c ">" DEGS.Trinity.fasta.transdecoder.blastpHits.pep`) usando uniprot
  * 581 unicos

Entonces:

```bash
# renombramos las cabeceras
pep=DEGS.Trinity.fasta.transdecoder.pfamHits.pep

awk '/^>/{print $1, $2; next}{print}' < $pep | sed -e 's/\::/ /' | awk '/^>/{print $1, $2; next}{print}' > DEGS.Trinity.fasta.pep

# read number of unique

grep ">" DEGS.Trinity.fasta.pep| sort | uniq | wc -l

```
>  Solo deberia haber un ORF por secuencia , sin embargo vemos que aun hay repeticiones

Preparamos el archivo de entrada para string

Debido a que tenemos mas de un marco de lectura para algunos transcritoos, re-etiquetamos las cabeceras de los fasta con numeros para evitar redundancia en los nombres en el programa string

```bash
awk '/^>/{print ">ORF_Sequence_" ++i; next}{print}' < DEGS.Trinity.fasta.pep > DEGS.Trinity.fasta.orfSeq.pep

```

### SuperTranscript

Dejemos de trabajar con isoformas y usemos genes, hay varias formas de trabajar:

* Usar la isoforma mas abundante (`$TRINITY_HOME/util/filter_low_expr_transcripts.pl`) como la representativa
* Usar la isoforma mas larga como la representativa
* Elaborar un scafold de todas las isoformas (SuperTranscript)

```bash
#Generate Trinity SuperTranscripts like so:
TRINITY_HOME=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/
$TRINITY_HOME/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta Trinity.fasta
```

Repetimos los primeros pasos de este documento.

```bash
awk '{print $1}' genes_DEGs.list | xargs -I {} grep -A1 '{}' trinity_genes.fasta > DEGS.trinity_genes.fasta
```

**Resultados**

* 135, 474 Genes (a lo largo del ensamble)
* DEGs (genes/isoformas):
  * **670** / 766 (`cut -d'_' -f1-4  DEGS.list | sort | uniq | wc -l`)
  * 603 DEGS.trinity_genes.fasta.transdecoder.pfamHits.pep 
  * 1, 889 longest ORFs

```bash
#
srun $DECODER/TransDecoder.LongOrfs -G universal -t DEGS.trinity_genes.fasta > DEGS_transdecoder.log &
#
srun $HMMSCAN/hmmscan --cpu 24 --domtblout DEGsPFAM.out ../../../Pfam-A.hmm longest_orfs.pep
#
srun blastp -query longest_orfs.pep      -db ../../../uniprot_sprot.pep  -max_target_seqs 1     -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6 &

# 
$DECODER/TransDecoder.Predict -t DEGS.trinity_genes.fasta --retain_pfam_hits DEGS.trinity_genes.fasta.transdecoder_dir/DEGsPFAM.out --single_best_orf

```

#plot like this
https://www.researchgate.net/figure/Length-histograms-of-the-Tritonia-Trinity-assembled-transcripts-and_fig9_271134380



## STRING

Evaluemos los marcos de lectura en string (https://string-db.org/)

- DEGs entran a transdecoder 
- Y despues la red de los DEGs con marco de lectura abierto
- tabla de DEGs en string y score >= 0.7 

comparar similitud de secuencia que interactuan en la red (String) con bases de datos de organismos modelo (humano) 

## igraph en R

```R
dir <- getwd()
setwd(dir)
k <- read.table("string_interaction_zero_score.tsv", header=FALSE)
colnames(k) <- c("node1", "node2", 
        "node1_string_internal_id",
        "node2_string_internal_id",
        "node1_external_id",       
        "node2_external_id",
        "neighborhood_on_chromosome",
        "gene_fusion",
        "phylogenetic_cooccurrence",       
        "homology",        
        "coexpression",      
        "experimentally_determined_interaction", 
        "database_annotated",     
        "automated_textmining",    
        "combined_score")

links <- k[,c("node1", "node2","combined_score")]


#::: prepare list to compare versus metadata
list <- c(as.vector(links$node2), as.vector(links$node1))
list <-list[order(list)]
list <- list[!duplicated(list)]

```

And load data:

```R
library(igraph)
net <- graph_from_data_frame(d = links, # nodes interaction
                             # vertices = nodes # attribudes per node
                             directed=TRUE)
```

The description of an [igraph](http://igraph.org/) object starts with four letters:

1. D or U, for a directed or undirected graph
2. N for a named graph (where nodes have a `name` attribute)
3. W for a weighted graph (where edges have a `weight` attribute)
4. B for a bipartite (two-mode) graph (where nodes have a `type` attribute)

The two numbers that follow (17 49) refer to the number of nodes and edges in the graph. The description also lists node & edge attributes, for example:

- `(g/c)` - graph-level character attribute
- `(v/c)` - vertex-level character attribute
- `(e/n)` - edge-level numeric attribute

```r
# Access vertices and edges:

E(net) # The edges of the object
V(net) # The vertices of the object

# You can also examine the network matrix directly:

net[]
net[1,] 

```

Lets visualize some properties of the net:

```R
plot(degree_distribution(net),
    #log="xy",
    col="sienna4",
    xlab = "Number of connections",
    ylab = "Density")
# or
boxplot(degree(net), col="snow2", ylab = 'Degree') # N of connections

# or
hist(degree(net, mode="all"), col = "tomato", 
    xlab = "Degree", breaks=20)


```

### Improve visualization

Notice that our network plot is still not too helpful. We can identify the type and size of nodes, but cannot see much about the structure since the links we’re examining are so dense. One way to approach this is to see if we can sparsify the network, keeping only the most important ties and discarding the rest.

```r
# And statiscs
hist(links$combined_score, breaks=20, col = "gray80", xlab = 'Score')
mean(links$combined_score)
sd(links$combined_score)
```

There are more sophisticated ways to extract the key edges, but for the purposes of this exercise we’ll only keep ones that have weight higher than the mean for the network. In igraph, we can delete edges using `delete_edges(net, edges)`:

```r
# 
# cut.off <- mean(links$combined_score)
# mean = 0.8888874
# (high confidence) - ACTUALLY REMOVED FROM STRING 
cut.off <- 0.7 # 
net <- delete_edges(net, E(net)[combined_score < cut.off])
```

>  Community detection (by optimizing modularity over partitions):

```r
# clp <- igraph::cluster_optimal(net)
#class(clp)
# At optimal_modularity.c:85 : GLPK is not available, Unimplemented function call
```

We could also use `simplify` to combine multiple edges by summing their weights with a command like `simplify(net, edge.attr.comb=list(Weight="sum","ignore"))`. Note, however, that this would also combine multiple edge types (in our data: “hyperlinks” and “mentions”).

```R
# coloring nodes
library("RColorBrewer")

cols <- brewer.pal(4, "Dark2")

V(net)$color <- ifelse(degree(net) == 0, cols[1], 
                       ifelse(degree(net) == 1, cols[2], 
                          ifelse(degree(net) > 10 & degree(net) < 20, 
                              cols[3], cols[4])))

V(net)$color <- ifelse(degree(net) == 1, cols[1], 
                       ifelse(degree(net) > 1 & degree(net) < 20, 
                              cols[2], cols[3]))

# Simplifing (if plot doesn't look very good)
net <- igraph::simplify(net, remove.multiple = F, remove.loops = T) 

plot(net, layout=layout_with_kk,
     edge.arrow.size=.2, # reduce arrow size
     vertex.size=5,
    vertex.label=NA, # remove labels
    #edge.curved=.1
    ) 
# using curved edges will allow you to see multiple links
# between two nodes (e.g. links going in either direction, or multiplex links) (timeDemand)

legend(x=-1.5, y=-1.1, c("None Connection",
                         "Single Connection",
                         "Multiple Connections",
                        "High Degree of Connections (>20)"), pch=21,
       col="#777777", pt.bg=cols, pt.cex=2, cex=.8, bty="n", ncol=1)
```

And generate table

```R
save <- data.frame(ifelse(degree(net) == 1, 'Single', 
                       ifelse(degree(net) > 1 & degree(net) < 20, 
                              'Multiple', 'High_Degree')), degree(net))


save <- data.frame(ifelse(degree(net) == 0, 'none', 
                       ifelse(degree(net) == 1, 'Single', 
                          ifelse(degree(net) > 1 & degree(net) < 20, 
                              'Multiple', 'High_Degree'))), degree(net))

names(save) <- c('Label','Degree_of_Connection')

write.csv(save, file = "Degree_of_Connections.csv", 
          row.names = TRUE, quote = FALSE)
```

or improve annotation as follow:

```R
string_annot <- read.csv('string_network_coordinates _zero_score.txt', header=FALSE, sep='\t')

names(string_annot) <- c('node',	'x_position',
                         'y_position',	'color',
                         'annotation')
```

Then,

```R
library(tidyverse)

save %>%
    as.tibble(rownames = 'node') %>%
    group_by(node) %>%
    inner_join(string_annot, by = "node") %>%
    select(node, Label, Degree_of_Connection, annotation) %>%
	write.csv(file = "Degree_of_Connections_zero_score.csv", 
          row.names = FALSE, quote = FALSE)
```

## Add Size of node

Let’s say we want to color our network nodes based on the degree attributes (more links -> larger node) We will also change the width of the edges based on their weight.

```r
# Compute node degrees (#links) and use that to set node size:
deg <- degree(net, mode="all")
V(net)$size <- deg*3

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label <- NA

# Set edge width based on combined_score:
E(net)$width <- E(net)$combined_score/6

#change arrow size and edge color:
E(net)$arrow.size <- .2
E(net)$edge.color <- "gray80"

# We can even set the network layout:
graph_attr(net, "layout") <- layout_with_lgl

plot(net,
     vertex.size=5)
```



# Coloring

```R
plot(x=1:5, y=rep(5,5), pch=19, cex=12, col=rgb(.25, .5, .3, alpha=.5), xlim=c(0,6))
#
palf <- colorRampPalette(c("gray70", "dark red", "orange"))  
plot(x=10:1, y=1:10, pch=19, cex=5, col=palf(10))

# Using RColorBrewer palettes in plots:
library("RColorBrewer")
pal3 <- brewer.pal(8, "Blues")
plot(x=10:1, y=10:1, pch=19, cex=6, col=pal3)
```

## Ortologia

Comparamos la ortologia de aquellos genes que estan dentro y fuera de la red (**aun no me que clara la idea de esta comparación**)

- usando herramienta proteinortho <- revisar apuntes
- http://pl.postech.ac.kr/QuickParanoid/....

Reference : 

https://github.com/TransDecoder/TransDecoder/wiki



