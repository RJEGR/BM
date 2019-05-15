## Asignacion de mock-sanger con bases de datos

Preparemos las bases de datos, midori, coarbitrator (genBank) y BOLD, para probar con los datos de la mock 68 con el mejor el valor omega que se obtuvieron el numero de ASVs==ntaxa(morfológico) `OMEGA_A_1e_120_ASVs.fasta`; 

1. Standar BOLD

```bash
#!/bin/bash
# bold
db=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/BOLD/BOLD_public.ALL.* 

ln -s $db .

sbatch rdp_assign.sh OMEGA_A_1e_120_ASVs.fasta 99 BOLD_public.ALL.tax

# for i in $db; do unlink $i; done
```



2. Midori

```bash
#!/bin/bash
# midori
db_midori=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/midori_unique_DB-0.* 

ln -s $db_midori .

sbatch rdp_assign.sh OMEGA_A_1e_120_ASVs.fasta 99 midori_unique_DB-0.tax

# for i in $db_midori; do unlink $i; done

```

3. Coarbitrator

```bash
#!/bin/bash
# coarbitrator

db_co=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/co-arbitrator/Coarbitrator_COI_nuc_curated.*

ln -s $db_co

sbatch rdp_assign.sh OMEGA_A_1e_120_ASVs.fasta 99 Coarbitrator_COI_nuc_curated.tax

# for i in $db_co; do unlink $i; done
```

extra. Only bold-species

```bash
#!/bin/bash
db_sp=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/BOLD/BOLD_public_species.*

ln -s $db_sp .

sbatch rdp_assign.sh OMEGA_A_1e_120_ASVs.fasta 99 BOLD_public_species.tax


```

