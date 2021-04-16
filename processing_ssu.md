

```bash
#!/bin/bash
#SBATCH --job-name=ssu
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --mail-user=rgomez@cicese.edu.mx
#SBATCH --ntasks-per-node=20
#SBATCH -t 6-00:00:00

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

#######################
# setting work_directory
# and exporting tools
########################
cd $SLURM_SUBMIT_DIR
export $SLURM_NPROCS
#mothur=/LUSTRE/bioinformatica_data/genomica_funcional/bin/mothur/mothur
export PATH="$PATH:/LUSTRE/bioinformatica_data/genomica_funcional/bin/ssu-align-0.1.1/bin"
export MANPATH="$MANPATH:/LUSTRE/bioinformatica_data/genomica_funcional/bin/ssu-align-0.1.1/share/man"
export SSUALIGNDIR="/LUSTRE/bioinformatica_data/genomica_funcional/bin/ssu-align-0.1.1/share/ssu-align-0.1.1"

echo 'Separating bacterial sequence due to HMM model'

echo " "

mkdir ssu_pick_dir

ssu-prep -f -x ./*.fasta ssu_pick_dir $SLURM_NPROCS

ssu_start=`date | cut -d" " -f4,2,3`

./ssu_pick_dir.ssu-align.sh

echo "ssu-align starts at $ssu_start and ends at $(date | cut -d" " -f4,2,3)"
echo 'Continuing mothur analysis (:'
echo "

# save accession ids
grep "^>" ./ssu_pick_dir/*bacteria.fa > ssu.bacteria.accnos
grep "^>" ./ssu_pick_dir/*eukarya.fa > ssu.eukarya.accnos

mv *accnos ./ssu_pick_dir
# draw

ssu-draw -a ./ssu_pick_dir/ssu_pick_dir.bacteria.stk

exit

```



...

```bash

# save accession ids
scrp -r rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/loberas/ssu_pick_dir

cd ssu_pick_dir

tail -n100 ssu_pick_dir.ssu-align.sum

```

```bash
#
# Summary statistics:
#
# model or       number  fraction        average   average
# category      of seqs  of total         length  coverage    nucleotides
# ------------  -------  --------  -------------  --------  -------------
  *input*          1921    1.0000         251.28    1.0000         482712
#
  archaea             0    0.0000              -         -              0
  bacteria         1909    0.9938         251.49    0.9977         480102
  eukarya             0    0.0000              -         -              0
#
  *all-models*     1909    0.9938         251.49    0.9977         480102
  *no-models*        12    0.0062         125.33    0.0000           1503
#
# Speed statistics:
#
# stage      num seqs  seq/sec  seq/sec/model    nucleotides    nt/sec
# ---------  --------  -------  -------------  -------------  --------
  search         1921  1921.000        640.333         482712  482712.0
  alignment      1909  1909.000       1909.000         480102  480102.0
```

Prepare data to EDA

```bash
grep 'ASV' ssu_pick_dir.bacteria.hitlist > ssu_pick_dir.bacteria.hitlist.txt
/Users/cigom/metagenomics/Loberas_MG/fastq/dada2_asv/ssu_pick_dir

```

