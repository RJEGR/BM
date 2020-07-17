```bash
conda activate qiime2-2020.6

fasta_in=dada_ASVs.fasta
count_table=dada_ASVs_count.table
tax_table=dada_ASVs.w2pr2_worms_API02.wang.taxonomy

# 1) Import Fasta
qiime tools import \
  --input-path $fasta_in \
  --output-path ${fasta_in%.fasta}.qza \
  --type 'FeatureData[Sequence]'

# 2) Import Tax
# We import these data into QIIME 2 Artifacts. Ex. the reference taxonomy file (w2pr2_worms_API02.tax) is a tab-separated (TSV) file without a header, we must specify HeaderlessTSVTaxonomyFormat as the source format since the default source format requires a header.

qiime tools import \
  --input-path $tax_table \
  --output-path ${tax_table%.taxonomy}.qza \
  --type 'FeatureData[Taxonomy]' \
  --input-format  "HeaderlessTSVTaxonomyFormat" #"HeaderlessTSVTaxonomyFormat"

# 3) feature count to  biom > qza
# http://biom-format.org/documentation/biom_conversion.html

biom convert \
	-i $count_table \
	-o ${count_table%.table}.biom \
	--to-hdf5 --table-type="OTU table"


qiime tools import \
  --input-path ${count_table%.table}.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ${count_table%.table}.qza
  


 # Test dataset
 qiime taxa barplot \
  --i-table ${count_table%.table}.qza \
  --i-taxonomy  ${tax_table%.taxonomy}.qza \
  --m-metadata-file metadata-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
 
 
 # exit 0
 
 # test
qiime feature-table summarize \
  --i-table ${count_table%.table}.qza \
  --o-visualization ${count_table%.table}.qzv

# https://forum.qiime2.org/t/issue-importing-tsv-taxonomy-file/1873
qiime feature-table summarize \
  --i-table ${tax_table%.taxonomy}.qza \
  --o-visualization ${tax_table%.taxonomy}.qzv
  
  
   qiime taxa barplot \
  --i-table ${count_table%.table}.qza \
  --i-taxonomy  ${tax_table%.taxonomy}.qza \
  --m-metadata-file metadata-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
  
  
```

