```bash
#!/bin/bash
### Directivas
#SBATCH -p cicese
#SBATCH --job-name=blastp
#SBATCH --ntasks-per-node=24
#SBATCH -N 2
#SBATCH -t 6-00:00:00


echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"




```





```bash
Last login: Thu Oct 21 09:18:35 on console
cigom@ciceses-MacBook-Pro ~ % ssh rgomez@omica
rgomez@omica's password:
Last login: Tue Oct 19 10:55:54 2021 from 158.97.225.72
[rgomez@omica genomica_funcional]$
[rgomez@omica genomica_funcional]$ htop
[rgomez@omica genomica_funcional]$ ls
Alison_DE  Clarissa      erick   Laura                 MEH     mmartinez        paoline_home-bak  rgomez         silvia_work_dir
anadely    claudia       fausto  Laura_file            MG_12S  mothur           pavel             rigoberto      Tripp
bin        common_files  FQS     linux_shell_tutorial  MG_18S  ojuarez          pkgs              run20_dinobac  vincent
Blanca     diana         giulia  LuisFigueroa          MG_28S  oyster_datasets  PULPO_Project     sandra
cei        Edgar         Hueman  manuelm               MG_COI  paoline          RAD2020           scripts
[rgomez@omica genomica_funcional]$ ^C
[rgomez@omica genomica_funcional]$ ^C
[rgomez@omica genomica_funcional]$ ^C
[rgomez@omica genomica_funcional]$ ls
Alison_DE  Clarissa      erick   Laura                 MEH     mmartinez        paoline_home-bak  rgomez         silvia_work_dir
anadely    claudia       fausto  Laura_file            MG_12S  mothur           pavel             rigoberto      Tripp
bin        common_files  FQS     linux_shell_tutorial  MG_18S  ojuarez          pkgs              run20_dinobac  vincent
Blanca     diana         giulia  LuisFigueroa          MG_28S  oyster_datasets  PULPO_Project     sandra
cei        Edgar         Hueman  manuelm               MG_COI  paoline          RAD2020           scripts
[rgomez@omica genomica_funcional]$ cd rgomez/
[rgomez@omica rgomez]$
[rgomez@omica rgomez]$ ls
blast2go_cli_20190805                          Genomes     miRNAs                               oyster-rawdata  trm-158917.log
cephalopod_full_assembly_annot                 loberas     oktopus_full_assembly                tar.sh          trm-158942.err
curso_miRNAs_cei                               md5sum.log  oyster_full_assembly                 test_2019_an    trm-158942.log
FLafarga_CICESE_Abalone_20161108−01479.tar.gz  mg          oyster_full_assembly_mirnas_minning  trm-158917.err  WGCNA
[rgomez@omica rgomez]$ less tar.sh
[rgomez@omica rgomez]$ less md5sum.log
[rgomez@omica rgomez]$ less md5sum.log
[rgomez@omica rgomez]$ cd ../^C
[rgomez@omica rgomez]$ ls cephalopod_full_assembly_annot/
ASSEMBLY_COMPLETE               transdecoder.log
blastp.sh                       trinity_genes.fasta
peerj-04-1763-s005.pep          trinity_genes.fasta.transdecoder.bed
peerj-04-1763-s005.pep.phr      trinity_genes.fasta.transdecoder.cds
peerj-04-1763-s005.pep.pin      trinity_genes.fasta.transdecoder_dir
peerj-04-1763-s005.pep.pog      trinity_genes.fasta.transdecoder.gff3
peerj-04-1763-s005.pep.psd      trinity_genes.fasta.transdecoder.pep
peerj-04-1763-s005.pep.psi      trinity_genes.fasta.transdecoder.pfam.pep
peerj-04-1763-s005.pep.psq      trinity_genes.fasta.transdecoder.pfam_vs_peerj-04-1763-s005_blastp.outfmt6
slurm-154149.out                trinity_genes.gtf
Transcriptoma_Referencia.fasta
[rgomez@omica rgomez]$ ls cephalopod_full_assembly_annot/blastp.sh
cephalopod_full_assembly_annot/blastp.sh
[rgomez@omica rgomez]$ less cephalopod_full_assembly_annot/blastp.sh
[rgomez@omica rgomez]$ ls
blast2go_cli_20190805                          Genomes     miRNAs                               oyster-rawdata  trm-158917.log
cephalopod_full_assembly_annot                 loberas     oktopus_full_assembly                tar.sh          trm-158942.err
curso_miRNAs_cei                               md5sum.log  oyster_full_assembly                 test_2019_an    trm-158942.log
FLafarga_CICESE_Abalone_20161108−01479.tar.gz  mg          oyster_full_assembly_mirnas_minning  trm-158917.err  WGCNA
[rgomez@omica rgomez]$
[rgomez@omica rgomez]$ ls
blast2go_cli_20190805                          Genomes     miRNAs                               oyster-rawdata  trm-158917.log
cephalopod_full_assembly_annot                 loberas     oktopus_full_assembly                tar.sh          trm-158942.err
curso_miRNAs_cei                               md5sum.log  oyster_full_assembly                 test_2019_an    trm-158942.log
FLafarga_CICESE_Abalone_20161108−01479.tar.gz  mg          oyster_full_assembly_mirnas_minning  trm-158917.err  WGCNA
[rgomez@omica rgomez]$ ls oyster_full_assembly
annotation       diffExp  gosemSim        samples.file      slurm-153658.log  trimmomatic.sh  trinityStats.log
assembly.sh      ExN50    multiqc         samples.txt       supertranscript   Trinity.fasta
Diff_Exon_Usage  fastqc   quantification  slurm-153658.err  transrate         trinity.full
[rgomez@omica rgomez]$ ls oyster_full_assembly/annotation/
build.log             Pfam-A.hmm      Pfam-A.hmm.h3p    superTranscript.annot    uniprot_sprot.pep
full.assembly.anont   Pfam-A.hmm.h3f  rebuild.sh        Trinotate.sqlite         uniprot_sprot.pep.phr
good.assembled.annot  Pfam-A.hmm.h3i  slurm-156473.out  Trinotate.TaxonomyIndex  uniprot_sprot.pep.pin
hmmpress.log          Pfam-A.hmm.h3m  slurm-156474.out  Trinotate.UniprotIndex   uniprot_sprot.pep.psq
[rgomez@omica rgomez]$ less oyster_full_assembly/annotation/rebuild.sh
[rgomez@omica rgomez]$
#!/bin/sh

## Directivas
#SBATCH --job-name=set-db
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=1

#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*tar.gz
#
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*tar.gz.md5
wget ftp://ftp.ncbi.nlm.nih.gov/b.ast/db/nr.*tar.gz.md5

./rgomez/oyster-rawdata/method_v2/ANNOTATE/NCBIdb/nr/download.sh (END)

```





```bash
#!/bin/sh

## Directivas
#SBATCH --job-name=set-db
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=1

cat urls | sh
./rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/download.sh (END)
```

