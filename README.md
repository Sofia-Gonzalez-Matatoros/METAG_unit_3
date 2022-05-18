# METAG_unit_3. Sofía González Matatoros
## Fecha: 20/05/2022
## Ejercicio
Repeat all the steps with a viral metagenome from a human saliva sample (virome.zip in Unit 3 of Moodle). Compare different de novo assemblies options (try _--meta--) or different kmer values. You must perform at least 3 different assemblies. Write a brief summary describing the bioinformatic pipeline you have followed (trimming, decontamination, improve in quality, number of reads remove in each step, etc.). Compare different de novo assemblies with QUAST and choose the best base on the obtained metrics.

NOTE: In the quality filtering step, modify the MINLEN argument considering the original read length. Consider that reads with a minimum of 50% of the average original size are ok for subsequent analyses.

NOTE: Importantly, you do not have a reference genome for a metagenome.

### 1. Preprocesado
#### 1.1. Descargamos el dataset
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

Para ello primeramente vemos la longitud media de nuestras lecturas.
```
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' virome_R1.fastq
```
> 298.939
```
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' virome_R2.fastq
```
> 299.249

La longitud media de los contings es 298 aproximadamente. Considerando que las lecturas que presenten la mitad de esta longitud tendrán una calidad aceptable, MINLEN se establecerá en 149.

Vamos a probar con tres combinaciones de parámetros, SLIDINGWINDOW:4:15 MINLEN:36 (la recomendada por el manueal de Trimmomatic), SLIDINGWINDOW:4:20 MINLEN:70 (la utilizada en la práctica en aula) y SLIDINGWINDOW:4:20 MINLEN:149 (adaptada a la longitud de nuestras lecturas)
```
trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qfa_paired.fq virome_R1_qfa_unpaired.fq virome_R2_qfa_paired.fq virome_R2_qfa_unpaired.fq SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qfb_paired.fq virome_R1_qfb_unpaired.fq virome_R2_qfb_paired.fq virome_R2_qfb_unpaired.fq SLIDINGWINDOW:4:20 MINLEN:70

trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qfc_paired.fq virome_R1_qfc_unpaired.fq virome_R2_qfc_paired.fq virome_R2_qfc_unpaired.fq SLIDINGWINDOW:4:20 MINLEN:149
```
#### 1.6. Comprobamos si la calidad ha mejorado y cuánto ha mejorado con cada uno de los parámetros
```
mkdir virome_QF_Quality
fastqc virome_R1_qfa_paired.fq -o virome_QF_Quality
fastqc virome_R2_qfa_paired.fq -o virome_QF_Quality
fastqc virome_R1_qfb_paired.fq -o virome_QF_Quality
fastqc virome_R2_qfb_paired.fq -o virome_QF_Quality
fastqc virome_R1_qfc_paired.fq -o virome_QF_Quality
fastqc virome_R2_qfc_paired.fq -o virome_QF_Quality
```
**Calidad R1 SLIDINGWINDOW:4:15 MINLEN:36**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r1a.png)

**Calidad R1 SLIDINGWINDOW:4:20 MINLEN:70**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r1b.png)

**Calidad R1 SLIDINGWINDOW:4:20 MINLEN:149**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/ef840fa874ca888726676b2d185319a72a4db241/fotos.1.6/1.6/r1c.png)

**Calidad R2 SLIDINGWINDOW:4:15 MINLEN:36**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r2a.png)

**Calidad R2 SLIDINGWINDOW:4:20 MINLEN:70**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/main/fotos.1.6/1.6/r2b.png)

**Calidad R2 SLIDINGWINDOW:4:20 MINLEN:149**
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/ef840fa874ca888726676b2d185319a72a4db241/fotos.1.6/1.6/r2c.png)

Se observa que las mejores calidades se consiguen con SLIDINGWINDOW:4:20 MINLEN:149

#### 1.7. Descontaminación: eliminamos las lecturas que alineen con los genomas humano y phiX174 usando Bowtie2
```
cd ~/Documents/unit_3_tarea
cp ~/Documents/unit_3/

mv ~/Downloads/Index.zip ~/Documents/unit_3_tarea
cd ~/Documents/unit_3_tarea
unzip Index.zip
cd Index
mv *.bt2 ../
cd ..
bowtie2 -x human-phix174 -1 virome_R1_qfc_paired.fq -2 virome_R2_qfc_paired.fq --un-conc virome_clean.fq -S tmp.sam
```
Nótese que el índice human-phiX174 está subido a moodle.

#### 1.8. Contamos cuántas lecturas nos quedan tras la descontaminación
```
wc -l *clean.* | awk '{print $1/4}'
```
> 56474
> 56474
> 112948

### 2. Ensamblado con SPAdes
Vamos a probar hasta tres ensamblados

```
spades.py -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_default
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_meta_33 -k33
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_meta_55 -k55
```
#### 2.1. Contamos el número de contigs y scaffolds

```
grep -c '>' ./virome_spades_default/*.fasta
```
> ./virome_spades_default/before_rr.fasta:3619
> ./virome_spades_default/contigs.fasta:3584
> ./virome_spades_default/scaffolds.fasta:3515

```
grep -c '>' ./virome_spades_meta_33/*.fasta
```
> ./virome_spades_meta_33/before_rr.fasta:6731
> ./virome_spades_meta_33/contigs.fasta:6136
> ./virome_spades_meta_33/first_pe_contigs.fasta:8811
> ./virome_spades_meta_33/scaffolds.fasta:5962

```
grep -c '>' ./virome_spades_meta_55/*.fasta
```
> ./virome_spades_meta_55/before_rr.fasta:5208
> ./virome_spades_meta_55/contigs.fasta:4967
> ./virome_spades_meta_55/first_pe_contigs.fasta:6086
> ./virome_spades_meta_55/scaffolds.fasta:4864


#### 2.2. Longitud y cobertura de los contings y scaffolds
```
grep '>' -m 5 ./virome_spades_default/*.fasta
grep '>' -m 5 ./virome_spades_meta_33/*.fasta
grep '>' -m 5 ./virome_spades_meta_55/*.fasta
```
| Default contings | meta_33 contings | meta_55 contings |
| ------------- | ------------- | ------------- |
| >NODE_1_length_87493_cov_10.414005 | >NODE_1_length_87399_cov_19.309983 | >NODE_1_length_87421_cov_17.225030 |
| >NODE_2_length_68667_cov_7.370528 | >NODE_2_length_66117_cov_10.316552 | >NODE_2_length_69098_cov_12.816694 |
| >NODE_3_length_66211_cov_4.952349 | >NODE_3_length_65493_cov_14.475588 | >NODE_3_length_66139_cov_9.040585 |
| >NODE_4_length_51161_cov_6.060803 | >NODE_4_length_49898_cov_12.550065 | >NODE_4_length_48968_cov_8.896469 |
| >NODE_5_length_48968_cov_5.306464 | >NODE_5_length_48968_cov_9.994421 | >NODE_5_length_43711_cov_11.021669 |


| Default scaffolds | meta_33 scaffolds | meta_55 scaffolds |
| ------------- | ------------- | ------------- |
| >NODE_1_length_87493_cov_10.414005 | >NODE_1_length_87399_cov_19.309983 | >NODE_1_length_87421_cov_17.225030 |
| >NODE_2_length_68667_cov_7.370528 | >NODE_2_length_66117_cov_10.316552 | >NODE_2_length_69098_cov_12.816694 |
| >NODE_3_length_66211_cov_4.952349 | >NODE_3_length_65493_cov_14.475588 | >NODE_3_length_66139_cov_9.040585 |
| >NODE_4_length_51161_cov_6.060803 | >NODE_4_length_49898_cov_12.550065 | >NODE_4_length_48968_cov_8.896469 |
| >NODE_5_length_48968_cov_5.306464 | >NODE_5_length_48968_cov_9.994421 | >NODE_5_length_43711_cov_11.021669 |

#### 2.3. Assembly stats
```
chmod a+x contigstats.py
python contigstats.py virome_spades_*/*.fasta
```
| sample | contigs | min | max | mean | n50 | bases | non_standard_bases |
| ------------- | ------------- | ------------- |------------- |------------- |------------- |------------- |------------- |
| virome_spades_default/contigs.fasta | 3584 | 128 | 87493 | 1005 | 1080 | 3602288 | 0 |
| virome_spades_default/scaffolds.fasta | 3515 | 128 | 87493 | 1025 | 1127 | 3603056 | 768 |
| virome_spades_meta_21/contigs.fasta | 4967 | 56 | 87421 | 735 | 704 | 3648784 | 0 |
| virome_spades_meta_21/scaffolds.fasta | 4864 | 56 | 87421 | 751 | 733 | 3650620 | 1836 |
| virome_spades_meta_33/contigs.fasta | 6136 | 34 | 87399 | 687 | 691 | 4212765 | 0 |
| virome_spades_meta_33/scaffolds.fasta | 5962 | 34 | 87399 | 707 | 727 | 4216704 | 3939 |
| virome_spades_meta_55/contigs.fasta | 4967 | 56 | 87421 | 735 | 704 | 3648784 | 0 |
| virome_spades_meta_55/scaffolds.fasta | 4864 | 56 | 87421 | 751 | 733 | 3650620 | 1836 |

### 3.Comparación de las estrategias de ensamblado con QUAST
En primer lugar hay que descargarse Quast en la máquina virtual porque la interfaz Web está dando problemas.
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
Vamos a crear una carpeta con todo lo necesaria para quast para que sea más fácil de ejecutar. Sin embargo, para evitar duplicar datos vamos a crear hipervínculos a los ficheros con los contingsy scaffolds.
```
cd /home/bgm/Documents/unit_3_tarea/
mkdir quast
ln -rs ./virome_spades_default/contigs.fasta ./quast/contigs_default.fasta 
ln -rs ./virome_spades_default/scaffolds.fasta ./quast/scaffolds_default.fasta 
ln -rs ./virome_spades_meta_33/contigs.fasta ./quast/contigs_meta_33.fasta 
ln -rs ./virome_spades_meta_33/scaffolds.fasta ./quast/scaffolds_meta_33.fasta 
ln -rs ./virome_spades_meta_55/contigs.fasta ./quast/contigs_meta_55.fasta 
ln -rs ./virome_spades_meta_55/scaffolds.fasta ./quast/scaffolds_meta_55.fasta 
```
Ejecutamos quast
```
quast ./quast/contigs_default.fasta ./quast/scaffolds_default.fasta ./quast/contigs_meta_33.fasta ./quast/scaffolds_meta_33.fasta ./quast/contigs_meta_55.fasta ./quast/scaffolds_meta_55.fasta 
```
![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/640fbaaae33a7e22834cfee4b863195ad255733a/fotos.3/cumulative.png)

![Screenshot](https://github.com/Sofia-Gonzalez-Matatoros/METAG_unit_3/blob/640fbaaae33a7e22834cfee4b863195ad255733a/fotos.3/stats.png)

#### 3.1. ¿Cuál es la mejor estrategia de ensamblado?
Compare different de novo assemblies with QUAST and choose the best based on the obtained metrics (smaller number of contigs, higher N50, smaller L50, longest total assembly length, etc.).
