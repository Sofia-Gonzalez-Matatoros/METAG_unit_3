# METAG_unit_3. Sofía González Matatoros
## Fecha:
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
> PONER LAS GRÁFICAS DEL HTML
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
> INSERTAR IMÁGENES
Se observa que las mejores son las obtenidas con MINLEN 70 (protocolo profe) en vez del que están en general en el manual
#### 1.7. Eliminamos las lecturas que alineen con los genomas humano y phiX174 (Bowtie2)
```
cd ~/Documents/unit_3_tarea
bowtie2 -x human-phix174 -1 virome_R1_qfb_paired.fq -2 virome_R2_qfb_paired.fq --un-conc virome_clean.fq -S tmp.sam
```
Nótese que el índice human-phiX174 está subido a moodle.

#### 1.8. Contamos cuántas lecturas nos quedan
```
wc -l *clean.* | awk '{print $1/4}'
```
### 2. Ensamblado con SPAdes
```
spades.py -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_isolate
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_isolate -k33
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_isolate -k21
```
#### 2.1. Contamos el número de contigs y scaffolds:
```
grep -c '>' ./virome_isolate/*.fasta
```
```
grep -c '>' ./virome_careful/*.fasta
```
```
grep '>' -m 5 ./virome_careful/*.fasta
grep '>' -m 5 ./virome_isolate/*.fasta
```
```
grep '>' ./virome_careful/scaffolds.fasta
grep '>' ./virome_isolate/scaffolds.fasta
```
```
python contigstats.py virome_*/*.fasta
```

### 3.Comparación de las estrategias de ensamblado con QUAST
