Using [SuperTranscripts](https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts), we can explore differential transcript usage (DTU). Differential transcript usage analysis is complementary to differential gene expression (DGE) and differential transcript expression (DTE) analysis. For details on how DTU, DGE, and DTE compare, see ["Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences" by Soneson, Love, and Robinson; F1000 2016](https://f1000research.com/articles/4-1521/v2). A figure from this paper shown below illustrates the differences between them.

In short, DTU always involves DTE but not always DGE; in the case of isoform switching, it's possible for the overall output of gene expression to remain the same but different isoforms will be expressed. As shown in [Davidson et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1284-1), by applying [DEX-Seq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) to Supertranscripts, we can explore differential transcript usage via different transcript segments showing up with statistically significant differences in read coverage in response to some condition or treatment.

In addition to Trinity, to run this you will need to have the [STAR](https://github.com/alexdobin/STAR/releases) or [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml) aligner installed, in addition to the [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) software that's part of the [Subread](https://academic.oup.com/nar/article/41/10/e108/1075719/The-Subread-aligner-fast-accurate-and-scalable) package. Just ensure that the 'featureCounts' utility is available via your PATH setting.

The mini-pipeline for performing DTU analysis involves aligning reads the SuperTranscripts, counting the reads that align to 'exonic' regions, and performing an exon-level differential expression analysis. These steps are encapsulated within a Trinity script that wraps aligners STAR or HISAT2 along with DEXseq, and assuming you have these tools installed, you can run it like so:

```bash
# Generate Trinity SuperTranscripts like so:

$TRINITY_HOME/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
       --trinity_fasta Trinity.fasta
# and this should generate two output files:

# trinity_genes.fasta   :supertranscripts in fasta format
# trinity_genes.gtf     :transcript structure annotation in gtf format
```

Then:

```bash

DEXEQ_P=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/Analysis/SuperTranscripts/DTU/

STAR_P=/LUSTRE/bioinformatica_data/genomica_funcional/bin/STAR/
export PATH=$STAR_P:$PATH

module load gcc-5.4.0

srun $DEXEQ_P/dexseq_wrapper.pl --genes_fasta trinity_genes.fasta --genes_gtf trinity_genes.gtf --samples_file samples.file --aligner STAR --CPU 24 &

#
# tab-delim samples file with:
# sample_name  replicate_name /path/to/left.fq.gz  /path/to/right.fq.gz
...
```

> /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster_full_assembly/Diff_Exon_Usage