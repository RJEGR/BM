Implementamos Transdecoder para predecir el marco de lectura:

```bash
DECODER=/LUSTRE/apps/bioinformatica/TransDecoder-3.0.1/

#1 
input_fasta='file.fasta'

srun $DECODER/TransDecoder.LongOrfs -G universal -t $input_fasta  2> transdecoder.log &

#2

srun $DECODER/TransDecoder.Predict -t $input_fasta --cpu 24 >> transdecoder.log & 

# transdecoder is finished.  See output files file.fasta.transdecoder.*
```

## Including homology searches as ORF retention criteria): 

Para maximizar aún más la sensibilidad para capturar los ORF que pueden tener un significado funcional, se puede escanear todos los ORF en busca de homología con proteínas conocidas y retener todos esos ORF. Esto se puede hacer de dos formas populares: una búsqueda BLAST en una base de datos de proteínas conocidas y la búsqueda de PFAM para identificar dominios de proteínas comunes. En el contexto de TransDecoder, esto se hace de la siguiente manera:

Retener el mejor ORF por transcrito (`--single_best_orf`) basado en los resultados de dominios conservados evaluados con HMMER y la base de datos PFam (`retain_pfam_hits`) o blastp (`--retain_blastp_hits`) y la base de peptidos _uniprot_sprot_

```bash
cd trandecoder_dir

# Using PFAM Domains
# HMMSCAN || Pfam-A.hmm
HMMSCAN=/LUSTRE/bioinformatica_data/genomica_funcional/bin/hmmer-3.1b2-linux-intel-x86_64/binaries

# 1.
srun $HMMSCAN/hmmscan --cpu 24 --domtblout PFAM.out Pfam-A.hmm longest_orfs.pep &

# PFAM.out is the name of the outfile
# Pfam-A.hmm is a database 
# longest_orfs.pep is the predicted ORFs with transdecoder

# BLASTP || uniprot_sprot
# using uniprot_sprot.pep as database
srun blastp -query longest_orfs.pep  \
    -db uniprot_sprot.pep  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6 &

# retain_pfam_hits and single_best_orf:
srun 
$DECODER/TransDecoder.Predict -t $input_fasta --retain_pfam_hits *.transdecoder_dir/PFAM.out --single_best_orf

# mv *.transdecoder.pep *.transdecoder.pfamHits.pep

# or --retain_blastp_hits and and single_best_orf

$DECODER/TransDecoder.Predict -t $input_fasta --retain_blastp_hits *.transdecoder_dir/blastp.outfmt6 --single_best_orf

# mv *.transdecoder.pep *transdecoder.blastpHits.pep
```

