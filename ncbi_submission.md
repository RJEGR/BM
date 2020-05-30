**INFORMACION** 

0. Ingresamos a la plataforma de ncbi submission

A. Creamos primero el codigo de [**BioSample**](https://submit.ncbi.nlm.nih.gov/subs/biosample/SUB7386605/submitter) :

1. submitter, 
2. general info, 
3. sample type (usually invertebrate template is selected), 
4. attributes <-- be carefully here! **cada archivo que se cargara en la parte de DATOS-SRA debe ser enlistado en estos metadatos.**
5. Una lista de codigos con formato SAMNXXXXXXX se generara para cada metadato registrado.

B. Creamos el **BioProject** que y asociamos la informacion del Biosamples (Ej. SAMN12213361, etc)

1. https://submit.ncbi.nlm.nih.gov/subs/bioproject
2. Seleccionas el project type, sample scope, **release date** and further

**DATOS**

C. Cargamos bibliotecas crudas en el **SRA,** 

1. asociamos el id de nuestro *BioProject* (Ej. PRJNA552592)
2. Generamos los metadatos para cada biblioteca y asociamos a su id de *BioSamples* (ie. SAMNXXXXXXX)

\4. Cargamos el ensamble en **TSA**

1. Incluimos los ids del SRA (**SRRXXXXXX)**
2. Cargamos las bibliotecas en formato crudo usando el explorador Web via HTTP o Aspera Connect plugin
3. O (si mas de 300 archivos o 10 GB) via ascp:
4. ascp -i aspera.openssh -QT -l 100m -k1 -d ./reads/ [subasp@upload.ncbi.nlm.nih.gov](mailto:subasp@upload.ncbi.nlm.nih.gov):uploads/cgalindo_cicese.mx_8mb0LWvC

4.1. o via, ftp

1. **Navigate to the source folder** where the files for submission are;
2. **Establish an FTP connection** using the credentials below:Address:
3.  open ftp-private.ncbi.nlm.nih.gov Username: subftp Password: w4pYB9VQ
4. **Navigate to your account folder**: cd uploads/cgalindo_cicese.mx_zKLGzRQw
5. **Create a subfolder (required!)** with a meaningful name:mkdir new_folder
6. **Navigate to the target folder** you just created:cd new_folder
7. **Copy your files into the target folder**: put file_name



## Install aspera

```bash
# Download the aspera software

wget https://download.asperasoft.com/download/sw/connect/3.9.9/ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.tar.gz
 
# Decompress
tar -xvf ibm-aspera-connect-*-linux-g2.12-64.tar.gz
 
# Install it
sh ibm-aspera-connect-*-linux-g2.12-64.sh
 
# Setup software license
mkdir -p ~/.ssh; ln -s ~/.aspera/connect/etc/asperaweb_id_dsa.openssh ~/.ssh/
 
# Setup path
export PATH=~/.aspera/connect/bin:$PATH
 
# To make sure the path is automatically available once you login laster on, add the command to ~/.bashrc
echo export PATH=~/.aspera/connect/bin:\$PATH >> $HOME/.bashrc

 
# Make sure the ascp command is available now
which ascp
```

After make your Bioproject > Biosamples > SRA submission you can figure out this part in the step of files submission

```bash
# wating for firewall permission
ascp -i $PWD/aspera_coi.openssh -QT -l 100m -k1 -d $PWD/raw subasp@upload.ncbi.nlm.nih.gov:uploads/cgalindo_cicese.mx_8mb0LWvC/coi_folder

#  -T                              Disable encryption
# -l MAX-RATE                     Max transfer rate
# -i PRIVATE-KEY-FILE             Private-key file name (id_rsa)

# for 18S

ascp -i $PWD/aspera_18S.openssh -QT -l100m -k1 -d $PWD/raw_18S subasp@upload.ncbi.nlm.nih.gov:uploads/cgalindo_cicese.mx_8mb0LWvC

```

or via ft

```bash
# works!!!!!! but. colappse after all


ftp -i

open ftp-private.ncbi.nlm.nih.gov

#Username: subftp
#Password: w4pYB9VQ

cd uploads/cgalindo_cicese.mx_zKLGzRQw/coi_folder


prompt # switch Interactive mode off.

mput *.* # broke

put 012-X05-G44-COI-AMB_S64_L001_R1_001.fastq.gz
put 015-X05-B15-COI-AMB_S14_L001_R1_001.fastq.gz
put 015-X05-F38-COI-AMB_S25_L001_R1_001.fastq.gz
```

```bash
#!/bin/sh
#SBATCH --job-name=sra_submit
#SBATCH -N 1
#SBATCH --mem=20GB
#SBATCH --ntasks-per-node=12

SERVER='ftp-private.ncbi.nlm.nih.gov'
USER='subftp'
PASSW='w4pYB9VQ'


INPUTDIR=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/procesamiento_ecologico/diversity_by_cruises/submission_libs/raw

LOCALBASE=uploads/cgalindo_cicese.mx_zKLGzRQw	

# not working

for file in $(ls *gz); 
do; 
ftp -i -n $SERVER << END_SCRIPT
user $USER $PASSW
lcd $INPUTDIR
cd $LOCALBASE
put $FILE
bye
END_SCRIPT

done



```



 ## ribosomal RNA (rRNA), metazoan COX1 Submission

Details:

**Metazoan (multicellular animal) Mitochondrial COX1 submissions** must meet the following requirements:

- All sequences are from metazoan (multicellular animal) organisms.
- All sequences in the **FASTA** file contain only mitochondrial COX1 sequence. Flanking sequence should not be included.
- The following information must be provided regarding the organism: isolate or specimen-voucher. Mitochondrial genetic code if organism is not in the NCBI taxonomy database.

**Eukaryotic rRNA and rRNA-ITS submissions** must meet the following requirements:

- All sequences are eukaryotic
- All sequences in the FASTA file contain sequences from one of the following types: nuclear small or large subunit ribosomal RNA, nuclear internal transcribed spacer 1 or 2, nuclear rRNA-ITS region, or mitochondrial or chloroplast small or large subunit ribosomal RNA.



API_KEY Submission

```
36d803882b420020206e11528782085d9207

NCBI_API_KEY=74e0e4bf2d8eaa3bb742c46316dbafe12909
echo "NCBI_API_KEY=74e0e4bf2d8eaa3bb742c46316dbafe12909" >> $HOME/.bash_profile
```

**Steps**

1. Create **[Bioproject](https://submit.ncbi.nlm.nih.gov/about/bioproject-biosample/)** and associated it to SRA 

- Public description 

```
Main:
Relationship between environmental conditions and zooplankton community structure in the  southern Gulf of Mexico

Description:

Zooplankton play a pivotal role in sustain the majority of marine ecosystems. The distribution patterns and diversity of zooplankton are key information for understanding the functioning of those ecosystems. Here, a  marine section of the Gulf of Mexico were assessed using multi-species targeted locus approach in order to understand the zooplankton community dynamics
```



2. SRA submission of raw paired-end libraries [here:](https://submit.ncbi.nlm.nih.gov/subs/sra/)

- Important step here is set the release inmediatly or in a specific date upon publication

2.  

3. 



**Glossary**

https://www.ncbi.nlm.nih.gov/books/NBK54364/

ref

https://wiki.rc.hms.harvard.edu/display/O2/Aspera+to+download+NCBI+SRA+data

https://www.ncbi.nlm.nih.gov/books/NBK242625/

https://downloads.asperasoft.com/en/downloads/8?list

https://downloads.asperasoft.com/



https://github.com/CandiceChuDVM/RNA-Seq/wiki/Tutorial:-How-to-upload-your-data-to-the-Sequence-Read-Archive-(SRA)%3F