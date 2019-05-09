Los marcos de lectura predichos asi como el genoma del ostion, son el conjunto de datos que implementaremos para usar como referencia para encontrar sitios de union a miRNA (transcripto) o potenciales precursores de miRNA (en el genoma).

Procesamos el trancriptoma de todos los tejidos de ostión y predecimos marcos de lectura abiertos con transdecoder, retenemos ORFs que tengan hit con la base de datos pfam

```bash
# /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/miRNAs

# paths
HMMSCAN=/LUSTRE/bioinformatica_data/genomica_funcional/bin/hmmer-3.1b2-linux-intel-x86_64/binaries
DECODER=/LUSTRE/apps/bioinformatica/TransDecoder-3.0.1/
# 
cd trinity_genes.fasta.transdecoder_dir
# 1
srun $HMMSCAN/hmmscan --cpu 24 --domtblout trinity_genes.fasta.transdecoder_dir/longest_orfs_PFAM.out Pfam-A.hmm trinity_genes.fasta.transdecoder_dir/longest_orfs.pep &

# 2
srun $DECODER/TransDecoder.Predict -t trinity_genes.fasta --retain_pfam_hits trinity_genes.fasta.transdecoder_dir/longest_orfs_PFAM.out --single_best_orf &


# rename the sequences in your reference fasta, just remove everything after the space
awk '/^>Gene/{gsub(/[::]/, " "); print ">"$2; next}{print}' trinity_genes.fasta.transdecoder.pfam.cds > trinity_genes.fasta.transdecoder.pfam.headers.cds

awk '/^>/{gsub(/ /, "_"); print ""$1; next}{print}' mature.fa > mature_headers.fa

```

Iniciamos el analisis con diversas herramientas para microRNAs:

## microtar

```bash
# precursor genome site:
sbatch microtar.sh GCA_002022765.4_C_virginica-3.0_genomic.fna hairpin_collapsed2dna.fa
# mature 3p UTR interaction:
sbatch microtar.sh trinity_genes.3p_UTR.fasta mature_collapsed2dna.fa
# and complete ORFs
sbatch microtar.sh trinity_genes.fasta.transdecoder.cds mature_collapsed2dna.fa

# microtar.sh
query=$1
reference=$2
filename=$(echo ${query%.*}_vs_${reference%.*}.microtar.tsv)
outfile=$(echo ${filename%.*}.log)
./microtar -t $query -q $reference -f $filename 2> $outfile
```

## mirDeep2

Vamos a probar `mirdeep2` y comparar con `microtar`(resultados previamente obtenidos) . Despues de la instalacion de mirdeep2 () ejecutamos ;

>  A diferencia de microtar mirdeep2 (mapper.pl) es mas rapido , logramos anotar el genoma de c.virginica y el meta-transcriptoma-genes.cds desde el portatil en un tiempo de ejecucion < 3 minutos. Mientras que microtar , la ejecucion de meta-transcriptoma-genes.cds contra la base de datos de miRs maduros (mature.fa) ha demorado mas de 1 dia.

```bash
# 1.
# remove_white_space_in_id.pl first, then:
bowtie-build trinity_genes.fasta.transdecoder.pfam.headers.cds trinity_genes.fasta.transdecoder.pfam.headers.cds

# collapse sequence (it include the reformat identifier tag) option -m

awk '/^>/{print $1; next}{print}' mature.fa > mature.headers.fa
awk '/^>/{print $1; next}{print}' hairpin.fa > hairpin.headers.fa

# source ~/.bashrc
# 2. then map
file=hairpin.headers.fa
ref=GCA_002022765.4_C_virginica-3.0_genomic.fa

mapper.pl $file \
    -c -q -j -l 17 \
    -m -u -n -i \
    -p $ref \
    -s hairpin_collapsed2dna.fa \
    -t ${file%.fa}_vs_${ref%.cds}.arf \
    -v 2> mapping2.out

# 3. and also  mirdeep2 (error RNAfold)
miRDeep2.pl mature_collapsed2dna.fa \
    $ref \
    ${file%.fa}_vs_${ref%.cds}.arf  \
    none \
    none \
    none 2> report.log
    
# or
  
miRDeep2.pl mature_collapsed2dna.fa $ref ${file%.fa}_vs_${ref%.cds}.arf none mature.headers.fa hairpin_collapsed2dna.fa 2>report2.log

# maybe need use this format
# awk '/^>/{print $1; next}{print}' trinity_genes.fasta.transdecoder.cds > trinity_genes.fasta.transdecoder.headers.cds

```

## targetScan



Usamos ademas `targetscan`, obteniendo 3p_UTR regiones predichos con transdecoder

```bash
grep "three_prime_UTR" trinity_genes.fasta.transdecoder.gff3 > three_prime_UTR.bed
```

y en R recortamos algunas secuencias de tamano menor a 100 nucleotidos:

```R
UTR_3p_file <- list.files(path_mirna, pattern = "*bed", full.names = TRUE)
UTR_3p <- read.csv(UTR_3p_file, sep = '\t', header = FALSE)

UTR_3p$length <- UTR_3p[,5] - UTR_3p[,4]
UTR_3p_out <- UTR_3p[UTR_3p$length >= 100,]

hist(UTR_3p_out$length, breaks = 1000, 
                        main = paste0("the minLenght: ", min(UTR_3p_out$length), "\nthe maxLenght: ", max(UTR_3p_out$length)))

write.table(UTR_3p_out[, -c(ncol(UTR_3p_out))], file =  paste0(path_mirna, "/", "three_prime_UTR.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')

```

Finalmente `bedtools`

```bash
bedtools getfasta -fi ../full_assembly/trinity_genes.fasta -bed three_prime_UTR.bed  > trinity_genes.3p_UTR.fasta
```

> no es tan sencillo usar targetscan debido al formato que solicita para el archivo trinity_genes.3p_UTR.fasta

## miRanda

```bash
# 
miranda mature_collapsed2dna.fa trinity_genes.3p_UTR.fasta -sc 160 -en -25 -out mature_collapsed2dna_vs_trinity_genes.3p_UTR.out

# and using the ORFs - pfma based

miranda mature_collapsed2dna.fa test.trinity_genes.fasta.transdecoder.pfam.headers.cds -sc 160 -en -25 -out mature_collapsed2dna_vs_TEST_trinity_genes.fasta.transdecoder.pfam.headers.out

# then make table
# names: Seq1,Seq2,Tot Score,Tot Energy,Max Score,MaxEnergy,Strand,Len1,Len2,Positions
grep "^>>" mature_collapsed2dna_vs_trinity_genes.3p_UTR.out > mature_collapsed2dna_vs_trinity_genes.3p_UTR.R.txt

# results in cluster

sbatch miranda.sh mature_collapsed2dna.fa trinity_genes.fasta.transdecoder.pfam.headers.cds
```



## Pulpo como modelo cefalopodo



```bash
TRINITY_HOME=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/
DECODER=/LUSTRE/apps/bioinformatica/TransDecoder-3.0.1/
UTILS=/LUSTRE/apps/bioinformatica/trinityrnaseq/util/support_scripts
HMMSCAN=/LUSTRE/bioinformatica_data/genomica_funcional/bin/hmmer-3.1b2-linux-intel-x86_64/binaries
DECODER=/LUSTRE/apps/bioinformatica/TransDecoder-3.0.1/

#1 meta-assembly of oyster samples tissues

#2 transform to super-transcript each isoforms-gene group
srun $TRINITY_HOME/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta Transcriptoma_Referencia.fasta &

#3 predict ORF from isoforms-gene group with transdecoder
srun $DECODER/TransDecoder.LongOrfs -G universal -t trinity_genes.fasta > transdecoder.log &

#
cd trinity_genes.fasta.transdecoder_dir

# 3.1 homology with pfam db
srun $HMMSCAN/hmmscan --cpu 24 --domtblout longest_orfs_PFAM.out Pfam-A.hmm longest_orfs.pep &

# 4 get orf based on best pfam homoly
srun $DECODER/TransDecoder.Predict -t trinity_genes.fasta --retain_pfam_hits trinity_genes.fasta.transdecoder_dir/longest_orfs_PFAM.out --single_best_orf &
# sbatch hmmscan.sh Pfam-A.hmm longest_orfs.pep

# 5. rename the sequences in your reference fasta, just remove everything after the space
awk '/^>Gene/{gsub(/[::]/, " "); print ">"$2; next}{print}' < trinity_genes.fasta.transdecoder.pep > trinity_genes.fasta.transdecoder.pfam.pep

# 5 BLASTP
query=trinity_genes.fasta.transdecoder.pfam.pep
reference=peerj-04-1763-s005.pep
DBS=$PWD

$BLAST/blastp -query $query \
         -db $DBS/$reference -num_threads 24 \
         -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen" -evalue 1e-5 \
          > $filename
          
# /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/cephalopod_full_assembly

# load it in R, and plot:


```

### From Kenny et al 2015 (x)

[here](https://www.sciencedirect.com/science/article/pii/S1874778715300088?via%3Dihub#ec0005)

Genomes were compared to all known [metazoan](https://www.sciencedirect.com/topics/earth-and-planetary-sciences/metazoan) miRNA sequences, as downloaded from [miRbase](https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/mirbase) ([Griffiths-Jones et al., 2008](https://www.sciencedirect.com/science/article/pii/S1874778715300088?via%3Dihub#bb0075)) on the 6th of February 2013 using BLASTN with the following settings: − word_size 11 -reward 5 -penalty − 4 -gapopen 8 -gapextend 6. 

Putative miRNA sequences obtained by blast were checked to confirm that both arms of the putative miRNA were present, for general [homogeneity](https://www.sciencedirect.com/topics/earth-and-planetary-sciences/homogeneity) of 5′ (seed) sequence, for robust hairpin structures and for a lack of similarity to known protein, tRNA or rRNA sequences. These criteria are similar to those found in [Tarver et al. (2012)](https://www.sciencedirect.com/science/article/pii/S1874778715300088?via%3Dihub#bb0265) and [Quah et al. (2015)](https://www.sciencedirect.com/science/article/pii/S1874778715300088?via%3Dihub#bb0190), although as small [RNA libraries](https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/rna-libraries) were not sequenced as a part of our investigation, criteria related to processing and overhang listed in [Tarver et al. (2012)](https://www.sciencedirect.com/science/article/pii/S1874778715300088?via%3Dihub#bb0265) were not included in our process. We also performed BLASTN comparison of identified contigs to the NCBI nr database to exclude the possibility that [contamination](https://www.sciencedirect.com/topics/earth-and-planetary-sciences/contamination) with human DNA (or DNA from other species represented in this database) underlay the identification of candidate miR loci.

> Genomic precursor sites

```bash
# 1. make blastdb

query=hairpin_collapsed2dna.fa
reference=GCA_002022765.4_C_virginica-3.0_genomic.fna
DBS=$PWD

makeblastdb -dbtype nucl -in $reference -out $reference

$BLAST/blastn -query $query \
	 -word_size 11 -reward 5 -penalty -4 \
	 -gapopen 8 -gapextend 6 \
         -db $DBS/$reference -num_threads $SLURM_NPROCS \
         -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen" -evalue 1e-5 \
          > $filename
          
# run in omica
sbatch blastn.sh hairpin_collapsed2dna.fa GCA_002022765.4_C_virginica-3.0_genomic.fna

```

does not make positive results using the mature_collapsed2dna reads **but good results using  harping_collapsed2dna reads;**

Further, use bowtie better or mapper.pl from mirdeep2 (previous used) `mature.headers_vs_GCA_002022765.4_C_virginica-3.0_genomic.fna.arf` results.

---

Usando bowtie:

```bash
# remove_white_space_in_id.pl first, then:
query=mature_collapsed2dna.fa
reference=GCA_002022765.4_C_virginica-3.0_genomic.fna
filename=$(echo ${query%.*}_vs_${reference%.*})

bowtie-build $reference $reference


#$SLURM_NPROCS

bowtie -p 12 -f -n 1 -e 80 -l 18 -a -m 1 --best --strata $reference --al ${filename}.align --un ${filename}.unalign $query 2> ${filename}.log

# -f query input files are (multi-)FASTA .fa/.mfa
# -n max mismatches in seed (can be 0-3, default: -n 2)
# -e max sum of mismatch quals across alignment 
# -l seed length for -n (default: 28)
# -a report all alignments per read
# -m suppress all alignments if > <int> exist (def: no limit)
# best hits guaranteed best stratum; ties broken by quality
# strata hits in sub-optimal strata aren't reported (requires --best)
# --al write aligned reads/pairs to a file 
# --un write unaligned reads/pairs to file
```

### from Wheeler B. et al, 2009, The deep evolution of metazoan microRNAs [metazoa- mollusc - halliotis - abalone]

Known miRNAs were annotated by identifying homologous mature and star miRNA sequences in miRBase (Griffith-Jones et al. 2007) release 10.1. **Standalone BLAST (blastn, version 2.2.27) was used to generate a list of candidate identities**. This list was filtered using three criteria, evaluated on an **ungapped global alignment** of the read and the hit sequence beginning at the 5' end: 

(1) sequence length must match within 2 nt; 

(2) positions 2–7 of the seed sequence must be identical; and 

(3) the remainder of the alignment may contain no more than three mismatches. Sequence reads that matched a known miRNA or miRNA star sequence within the above criteria were annotated and removed from the data set.

```bash
query=hairpin_collapsed2dna.fa
reference=GCA_002022765.4_C_virginica-3.0_genomic.fna
DBS=$PWD

makeblastdb -dbtype nucl -in $reference -out $reference

$BLAST/blastn -query $query \
	 -word_size 11 -reward 5 -penalty -4 \
	 -gapopen 0 -gapextend 0 \
         -db $DBS/$reference -num_threads $SLURM_NPROCS \
         -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen" -evalue 1e-5 \
          > $filename
```



#### find novel miRNAs 

Known non-miRNA reads were identified by comparison with NCBI’s ‘‘nt’’ nucleotide database using Standalone MEGABLAST (version 2.2.17). Reads matching a known RNA molecule with percent identity 495% were removed from the data set, and the remaining sequence reads were then investigated for phylogenetic conservation. Reads from all species were combined, and those that ‘‘matched’’ a read with a higher frequency count were grouped. Matches were determined using the three criteria used to identify known miRNAs given above (similar length, seed sequence identity, and nonseed sequence similarity). 

> find conserved miRNAs from different metazoan

Reads conserved across multiple taxa were grouped, and groups were ranked by the frequency count of the most frequently occurring sequence. Reads not conserved across multiple taxa were divided by taxon and ranked by frequency count. This completed the automated analysis by miRMiner, resulting in a list of conserved reads across all taxa and lists of unique reads for each taxon