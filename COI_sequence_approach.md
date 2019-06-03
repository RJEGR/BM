> From Heller P. et al. 2018, Data Descriptor: A database of metazoan cytochrome c oxidase subunit I gene sequences derived from GenBank with CO-ARBitrator

The success of molecular approaches promped to concept of a 'barcode of animal life' - a gene whose sequence reliably uniquely identifies most animals. In 2002 Hebert proposed cytochrome c oxidase I (COI or cox1) as a standard for molecular barcoding of animals (Evolution of genes and taxa: a primer and Biological identification through DNA barcodes). COI is present in all animals; it has no introns; indels mutation are rare; and it has a high subtituton rate in the 3rd codon position, providing **nucleotide sequence diversity**.



> from Pentinsaari M. et al 2016

COI (a standardized 658 bp fragment) was proposed as a universal marker for species identification - to be used as a 'DNA barcode' tagging any taxon in the animal kingdom. Following this seminal idea, the number of partial COI gene sequences available in public data repositories has skyrocketed. DNA barcoding studies published to date treat this gene region as a mere identification tag - in exact accordance with the concept of a conveniently readable 'barcode'. Yet, the barcode fragent is located at the core of energy production within cells: **COI is one of th building blocks of the cytochrome C oxidase protein (COX)**.

The COX protein is a dimer compossed of two identical parts. These, in turn, consist of several amino acid chains (11 nuclear-encoded and three mitochindrial - encoded in mammals) as well as several metallic ligand: two iron atoms bound in heme group, three coppers, one zing and one magnesium ([Tsukihara et al 1996](https://science.sciencemag.org/content/272/5265/1136) and [Balsa E et al 2012](https://www.sciencedirect.com/science/article/pii/S1550413112002938?via%3Dihub)). COX is the last enzyme in the electron transport chain, reducing oxygen and pumping protons across the inner mitochondrial membrane. Thus, **changes in the amino acid sequences that modify the protein structure may affect energy metabolism**.

Amino acid subtitutions are rare especially in the cytochrome oxidase genes. These selective constraints on the aminoacid sequence are reflected at the DNA sequence leve: **the DNA barcode sequence cannot vary freely and its evolution is far from neutral**. 



**Nucleotide substitution rate in mammalian mitochondrial genomes**: 

> ref <https://www.ncbi.nlm.nih.gov/pubmed/10079281>

(1) High intragenomic variability in the evolutionary dynamic of mtDNA was found. The substitution rate is strongly dependent on the region considered, and slow- and fast-evolving regions can be identified. Nonsynonymous sites, the D-loop central domain, and tRNA and rRNA genes evolve much more slowly than synonymous sites and the two peripheral D-loop region domains. The synonymous rate is fairly uniform over the genome, whereas the rate of nonsynonymous sites depends on functional constraints and therefore differs considerably between genes. 

(2) The commonly accepted statement that mtDNA evolves more rapidly than nuclear DNA is valid only for some regions, thus it should be referred to specific mitochondrial components. In particular, nonsynonymous sites show comparable rates in mitochondrial and nuclear genes; synonymous sites and small rRNA evolve about 20 times more rapidly and tRNAs about 100 times more rapidly in mitochondria than in their nuclear counterpart. 

(3) A species-specific evolution is particularly evident in the D-loop region. As the divergence times of the organism pairs under consideration are known with sufficient accuracy, absolute nucleotide substitution rates are also provided.



### TEST

- non-synonymous change and synonymous change (dN/dS) sequences analysis example [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0002201) and tutorial [here](https://www.researchgate.net/publication/268231340_A_Beginners_Guide_to_Estimating_the_Non-synonymous_to_Synonymous_Rate_Ratio_of_all_Protein-Coding_Genes_in_a_Genome).

```R

.cran_packages <- c("dplyr")
.bioc_packages <- c("Biostrings", "seqinr")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F, version = "3.8")
}


# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)



# install.packages("seqinr", repos="http://R-Forge.R-project.org")

# Ks and Ka are, respectively, the number of substitutions per synonymous site and per non-synonymous site between two protein-coding genes. 
# They are also denoted as ds and dn in the literature. 
# The ratio of nonsynonymous (Ka) to synonymous (Ks) nucleotide substitution rates is an 
# indicator of selective pressures on genes. 
# A ratio significantly greater than 1 indicates positive selective pressure. 
# A ratio around 1 indicates either neutral evolution at the protein level or an averaging of sites under positive and negative selective pressures. 
# A ratio less than 1 indicates pressures to conserve protein sequence (i.e. purifying selection). 
# This function estimates the Ka and Ks values for a set of aligned sequences using the method published by Li (1993) and gives the associated variance matrix.


# path <- '/Users/cigom/metagenomics/COI/run012/mock_P68/mock_parameter_definition'
path <- '/Users/cigom/metagenomics/COI/run012/mock_P68/mock_parameter_definition/nucleotide_to_aa/'
# path <- '/Users/cigom/metagenomics/COI/run012/mock_P68/mock_parameter_definition/nucleotide_to_aa/sanger_to_aa'
setwd(path)

fasta.file <-  list.files(pattern = "pfam.pep") # 'ictio_coi_sanger114.fasta.transdecoder.pfam.pep'
fasta.file <- "OMEGA_A_1e_120_ASVs.fasta"

alphabetFrequency(readDNAStringSet(list.files(pattern = "pfam.pep")), baseOnly=TRUE, as.prob=TRUE)

# Construct phylogenetic tree
seqs <- readDNAStringSet(fasta.file)

# https://web.stanford.edu/class/bios221/labs/biostrings/lab_1_biostrings.html
alphabetFrequency(seqs, baseOnly=TRUE)

GCstaph = data.frame(ID=names(seqs), GC=rowSums(alphabetFrequency(seqs)[, c(2,3)]/width(seqs))*100)
window = 100
# compute the GC content in a sliding window (as a fraction) for a sequence no. 364
gc = rowSums(letterFrequencyInSlidingView(seqs[[69]], window, c("G", "C")))/window
plot(gc, type = 'l')

plot(1:length(gc), gc)
lines(lowess(x = 1:length(gc), y= gc, f = 0.10), col = 12, lwd = 2)



seqs.width <- width(seqs)
names_seqs <- names(seqs)
names(seqs) <- seqs # This propagates to the tip labels of the tree

alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA)
alignment <- subseq(alignment, start=3, end=302)
width(alignment)

# seqinr

alignment <- seqinr::as.alignment(nb = length(alignment), seq = as.data.frame(alignment)$x, com = NULL)
# subseq(bstring, start=3, end=3)



res <- seqinr::kaks(alignment)

if(any(!is.finite(res$ka))) stop("Non finite value returned for Ka")
if(any(!is.finite(res$ks))) stop("Non finite value returned for Ks")

plot(res$ks, res$ka)

# mmm no sure what is it!

heatmap(as.matrix(res$ka))


```



- AT content (see Pentinsaari M. et al 2016 nucleotide and a.a statistics and figure 1) - Nucleotide frequencies at each codon position. 
- implement nr search to download CO1_COI_COX1_COXI_GENE_Eukaryota_nr.fasta (teresita et al over 2.5 million COI seqs)