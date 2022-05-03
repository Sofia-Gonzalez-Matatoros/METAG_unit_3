# METAG_unit_3. Sofía González Matatoros
## Fecha: 20/05/2022
## Ejercicio
Repeat all the steps with a viral metagenome from a human saliva sample (Virome.zip in Unit 3 of Moodle). Compare different de novo assemblies options (try _--meta--) or different kmer values. You must perform at least 3 different assemblies. Write a brief summary describing the bioinformatic pipeline you have followed (trimming, decontamination, improve in quality, number of reads remove in each step, etc.). Compare different de novo assemblies with QUAST and choose the best base on the obtained metrics.

NOTE: In the quality filtering step, modify the MINLEN argument considering the original read length. Consider that reads with a minimum of 50% of the average original size are ok for subsequent analyses.

NOTE: Importantly, you do not have a reference genome for a metagenome.

### 1. Preprocesado
#### 1.1. Descargar el dataset
```
mkdir -p ~/Documents/unit_3_tarea
mv ~/Downloads/virome.zip ~/Documents/unit_3_tarea
cd ~/Documents/unit_3_tarea
unzip virome.zip
```
#### 1.2. Comprobamos la integridad de los archivos
```
md5sum virome_R1.fastq 
```
> cd891dfc865f01c8f3923edd55dacde5  virome_R1.fastq
```
md5sum virome_R2.fastq 
```
> 28f0d0ace9fb45af8b2595cbc78fa2bd  virome_R2.fastq
#### 1.3. Contamos el número de lecturas
```
wc -l virome_R1.fastq | awk '{print $1/4}'
```
> 100000
```
wc -l virome_R2.fastq | awk '{print $1/4}'
```
> 100000
#### 1.4. Comprobamos la calidad de las secuencias con fastqc
```
mkdir virome_Quality
fastqc virome_R1.fastq -o virome_Quality/
fastqc virome_R2.fastq -o virome_Quality/
```
**Calidad de R1**
![Screeshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.4/1.4/virome_R1.png)

**Calidad de R2**
![Screeshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.4/1.4/virome_R2.png)

#### 1.5. Recortamos los extremos de baja calidad con Trimmomatic
```
trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qfa_paired.fq virome_R1_qfa_unpaired.fq virome_R2_qfa_paired.fq virome_R2_qfa_unpaired.fq SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qfb_paired.fq virome_R1_qf_unpaired.fq virome_R2_qfb_paired.fq virome_R2_qfb_unpaired.fq SLIDINGWINDOW:4:20 MINLEN:70
```
#### 1.6. Comprobamos si la calidad ha mejorado y cuánto ha mejorado con cada uno de los parámetros
```
mkdir virome_QF_Quality
fastqc virome_R1_qfa_paired.fq -o virome_QF_Quality
fastqc virome_R2_qfa_paired.fq -o virome_QF_Quality
fastqc virome_R1_qfb_paired.fq -o virome_QF_Quality
fastqc virome_R2_qfb_paired.fq -o virome_QF_Quality
```
**Calidad R1 SLIDINGWINDOW:4:15 MINLEN:36**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r1a.png)

**Calidad R1 SLIDINGWINDOW:4:20 MINLEN:70**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r1b.png)

**Calidad R2 SLIDINGWINDOW:4:15 MINLEN:36**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r2a.png)

**Calidad R2 SLIDINGWINDOW:4:20 MINLEN:70**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r2b.png)

Se observa que las lecturas con mejor calidad son las obtenidas con el protocolo de prácticas (SLIDINGWINDOW:4:20 MINLEN:70), en comparación con las del manual de Trimmomatic (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

#### 1.7. Eliminamos las lecturas que alineen con los genomas humano y phiX174 (Bowtie2)
```
cd ~/Documents/unit_3_tarea
cp ~/Documents/unit_3/

mv ~/Downloads/Index.zip ~/Documents/unit_3_tarea
cd ~/Documents/unit_3_tarea
unzip Index.zip
cd Index
mv *.bt2 ../
cd ..
bowtie2 -x human-phix174 -1 virome_R1_qfb_paired.fq -2 virome_R2_qfb_paired.fq --un-conc virome_clean.fq -S tmp.sam
```
Nótese que el índice human-phiX174 está subido a moodle.

#### 1.8. Contamos cuántas lecturas nos quedan
```
wc -l *clean.* | awk '{print $1/4}'
```
> 80761
> 
> 80761
> 
> 161522

### 2. Ensamblado con SPAdes
```
spades.py -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_default
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_meta_33 -k33
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_meta_21 -k21
```
#### 2.1. Contamos el número de contigs y scaffolds

```
grep -c '>' ./virome_spades_default/*.fasta
```
> ./virome_spades_default/before_rr.fasta:4026
> ./virome_spades_default/contigs.fasta:3973
> ./virome_spades_default/scaffolds.fasta:3886

```
grep -c '>' ./virome_spades_meta_33/*.fasta
```
> ./virome_spades_meta_33/before_rr.fasta:7554
> ./virome_spades_meta_33/contigs.fasta:6687
> ./virome_spades_meta_33/first_pe_contigs.fasta:10295
> ./virome_spades_meta_33/scaffolds.fasta:6425

```
grep -c '>' ./virome_spades_meta_21/*.fasta
```
> ./virome_spades_meta_21/before_rr.fasta:9002
> ./virome_spades_meta_21/contigs.fasta:7421
> ./virome_spades_meta_21/first_pe_contigs.fasta:14591
> ./virome_spades_meta_21/scaffolds.fasta:7142

#### 2.2. Longitud y cobertura de los contings y scaffolds
```
grep '>' -m 5 ./virome_spades_default/*.fasta
grep '>' -m 5 ./virome_spades_meta_33/*.fasta
grep '>' -m 5 ./virome_spades_meta_21/*.fasta
```
| Default contings | meta_33 contings | meta_21 contings |
| ------------- | ------------- | ------------- |
| >NODE_1_length_87493_cov_11.125083 | >NODE_1_length_87399_cov_21.743287 | >NODE_1_length_87379_cov_23.132524 |
| >NODE_2_length_70174_cov_7.587063 | >NODE_2_length_65568_cov_15.579110 | >NODE_2_length_59788_cov_14.221627 |
| >NODE_3_length_66299_cov_6.577480 | >NODE_3_length_50892_cov_16.263198 | >NODE_3_length_48915_cov_12.334622 |
| >NODE_4_length_51626_cov_4.664129 | >NODE_4_length_49862_cov_16.341187 | >NODE_4_length_35664_cov_80.708554 |
| >NODE_5_length_51161_cov_7.061743 | >NODE_5_length_48968_cov_11.565791 | >NODE_5_length_34551_cov_16.385172 |


| Default scaffolds | meta_33 scaffolds | meta_21 scaffolds |
| ------------- | ------------- | ------------- |
| >NODE_1_length_87493_cov_11.125083 | >NODE_1_length_87399_cov_21.743287 | >NODE_1_length_87379_cov_23.132524 |
| >NODE_2_length_70174_cov_7.587063 | >NODE_2_length_65568_cov_15.579110 | >NODE_2_length_66226_cov_16.899857 |
| >NODE_3_length_66299_cov_6.577480 | >NODE_3_length_50892_cov_16.263198 | >NODE_3_length_59788_cov_14.221627 |
| >NODE_4_length_51626_cov_4.664129 | >NODE_4_length_49862_cov_16.341187 | >NODE_4_length_48915_cov_12.334622 |
| >NODE_5_length_51161_cov_7.061743 | >NODE_5_length_48968_cov_11.565791 | >NODE_5_length_40347_cov_9.057878 |

#### 2.3. Assembly stats
```
chmod a+x contigstats.py
python contigstats.py virome_spades_*/*.fasta
```
| sample | contigs | min | max | mean | n50 | bases | non_standard_bases |
| ------------- | ------------- | ------------- |------------- |------------- |------------- |------------- |------------- |
virome_spades_default/before_rr.fasta | 4026 | 128 | 67569 | 1031 | 1172 | 4150192 | 0 |
virome_spades_default/contigs.fasta | 3973 | 128 | 87493 | 1044 | 1200 | 4146223 | 0 |
virome_spades_default/scaffolds.fasta | 3886 | 128 | 87493 | 1067 | 1284 | 4147575 | 1352 |
virome_spades_meta_21/before_rr.fasta | 9002 | 22 | 33760 | 544 | 611 | 4901025 | 0 |
virome_spades_meta_21/contigs.fasta | 7421 | 22 | 87379 | 660 | 752 | 4901404 | 0 |
virome_spades_meta_21/first_pe_contigs.fasta | 14591 | 22 | 87136 | 476 | 772 | 6944615 | 0 |
virome_spades_meta_21/scaffolds.fasta | 7142 | 22 | 87379 | 688 | 824 | 4910698 | 9313 |
virome_spades_meta_33/before_rr.fasta | 7554 | 34 | 87399 | 631 | 649 | 4766605 | 0 |
virome_spades_meta_33/contigs.fasta | 6687 | 34 | 87399 | 710 | 756 | 4747439 | 0 |
virome_spades_meta_33/first_pe_contigs.fasta | 10295 | 34 | 75432 | 649 | 938 | 6685342 | 0 |
virome_spades_meta_33/scaffolds.fasta | 6425 | 34 | 87399 | 740 | 817 | 4756939 | 9500 |


### 3.Comparación de las estrategias de ensamblado con QUAST
Nos descargamos Quast
```
cd /home/bgm/
mkdir software
cd software # This folder will store download programs
wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
tar -xzf quast-5.0.2.tar.gz
rm quast-5.0.2.tar.gz
cd quast-5.0.2
python setup.py install # for mac
cd /home/bgm/
mkdir bin
cd bin # This folder will store symbolic links to program executables
ln -s /home/bgm/software/quast-5.0.2/quast.py quast
# edit /home/bgm/.bashrc with _vim_ or _nano_ and add the following line
export PATH=/home/bgm/bin/:$PATH
# This way, all executable that you link in the bin folder will be in the system path
# Open a new terminal window/tab
cd /home/bgm/software/quast-5.0.2/quast_libs/site_packages/jsontemplate
# Always make a copy of the file you are going to modify
cp jsontemplate.py jsontemplate_original.py
sed s/cgi/html/g jsontemplate_original.py > jsontemplate.py
```
Creamos hipervínculos para no duplicar los datos
```
cd /home/bgm/Documents/unit_3_tarea/
mkdir quast
ln -rs ./virome_spades_default/contigs.fasta ./quast/contigs_default.fasta 
ln -rs ./virome_spades_default/scaffolds.fasta ./quast/scaffolds_default.fasta 
ln -rs ./virome_spades_meta_33/contigs.fasta ./quast/contigs_meta_33.fasta 
ln -rs ./virome_spades_meta_33/scaffolds.fasta ./quast/scaffolds_meta_33.fasta 
ln -rs ./virome_spades_meta_21/contigs.fasta ./quast/contigs_meta_21.fasta 
ln -rs ./virome_spades_meta_21/scaffolds.fasta ./quast/scaffolds_meta_21.fasta 
```
Ejecutamos quast
```
quast ./quast/contigs_default.fasta ./quast/scaffolds_default.fasta ./quast/contigs_meta_33.fasta ./quast/scaffolds_meta_33.fasta ./quast/contigs_meta_21.fasta ./quast/scaffolds_meta_21.fasta 
```
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/0aef989106a94929e402ffbf780898665a2ee3df/Captura%20de%20pantalla%20de%202022-05-03%2019-49-33.png)

![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/0aef989106a94929e402ffbf780898665a2ee3df/Captura%20de%20pantalla%20de%202022-05-03%2019-49-43.png)
