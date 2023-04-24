

# Introduction 

We received sequencing data from the illumina NovaSeq 6000 processing 42 samples. 
In the script, sockeye preparations we have downloaded the tools we will use for most of the pipeline, as well as prepared an index for genome alignment. 

Paired reads, number of samples per lane: 

## Aims

- Quality control of sequencing (raw reads)
- Align them to current assembly of human genome (GENCODE v38) using STAR
- Summarize reads using both RSEM and feature counts
- Perform QC of alignment using picard tools
- Export data for differential gene expression in R

# Preparation steps

## Installing necessary software/data

### Reference genome

GRCh38.p13 v 38 (20210521) from GENCODE: 

- __GRCh38.p13.genome.fa.gz__ includes reference chromosomes, scaffolds, assembly patches and haplotypes.
- __gencode.v38.annotation.gtf.gz__ contains annotation for reference genome 
- Download using Rsync, unpack and decompress

```sh
cd target_directory

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v38.annotation.gtf.gz

# unpack in background
gunzip *.gz &
```

### STAR v 2.7.8

Download data source in Github, unpack and unzip. Compile using `make`

```sh
mkdir ~/software
cd ~/software/

wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
tar -xzf 2.7.8a.tar.gz
cd STAR-2.7.8a/source 
make STAR
```

### RSEM v 1.3.3

Same steps: download, unpack and make locally 

```sh
wget https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz

tar -xzf RSEM*
cd <PATH to RSEM>/bin | make
```

### FastQC v 0.11.9

Download and unpack `fastqc` from the binary source. Run `make` inside unpacked directory to create executable

```sh
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 755 fastqc
```

### Salmon v1.5.2

```sh
cd ~/software
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz
tar xzvf salmon-1.5.2_linux_x86_64.tar.gz
# salmon can be run from the bin folder

```

### Picard tools v2.26

Requires java

```sh
java -version
cd ~/software
wget https://github.com/broadinstitute/picard/releases/download/2.26.9/picard.jar
java -jar /path/to/picard.jar -h 

# download a reference flat file for the genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz

```

### MultiQC v 1.11

Clone git repository and install `python::pip` 

```sh
#in cluster
module load python
module load git 

cd ~/software
git clone https://github.com/ewels/MultiQC.git
cd MultiQC
pip install .
```

### Load the RNAseq data

Data for reproducibility must be obtained from the NCBI biorepository - https://www.ncbi.nlm.nih.gov/bioproject/PRJNA844016

Rename files to a suitable name `CFBXXX_R1.fastq.gz` for your local jobs


# Quality control of sequencing output


```sh

# run from scratch directory, create output dir and establish path to fastqc
cd ~
FASTQC='/path/to/fastqx/tool'
[ -d ~/QC_ntm/fastqc_raw ] || mkdir ~/QC_ntm/fastqc_raw # create only if it does not exist

# loop for all raw sequences
for fastq in ~/NTM_fasta/CFB*
do 
echo running "$fastq"
$FASTQC/fastqc -o QC_ntm/fastqc_raw \
                -t 8 \
                $fastq

done
```

## FastQC results

* Uniform distribution of nucleotides except in the start of reads, _Expected for RNA_seq_ as random primers are not completely random
* Non-normal GC distribution over reads, can be caused by _duplicated/adaptor sequences_
* High percentage of overrepresented sequences, probably _hemoglobin subunits or adaptor sequences_ 

# Alignment/quantification 

* **RSEM** aligns to the transcriptome and summarizeds reads to transcript or gene level. 
* Same coverage as novel algorithms ( **salmon|kallisto** ) , but more time-consuming.
* More established and comparable to previous results. 
* STAR is an ultrafast aligner to a reference genome

## Strandness of reads

Using the `salmon` pseudoaligner, we infer the sense and strandness of reads. Only necessary to run a pilot subset of samples. 

Necessary for quantification and picard quality control after alignment.

__Library is type ISR.__

- I= Inward direction, the paired reads point towards the same point, inside the sequence
- S= Stranded, there is a stranded library preparation
- R= Means that read 1 comes from the reverse strand

```sh
# create variables with paths to GENCODE and salmon
salmon_dir='/path/to/fastqx/tool'
GENCODE='/path/to/reference/genome'
SCRATCH='/path/to/working/directory'
fasta='~/NTM_fasta'

# create transcriptome index
$salmon_dir/salmon index \
	-t $GENCODE/gencode.v38.transcripts.fa \
    -i $SCRATCH/gencode_v38_salmon_idx
							

# while loop, starts with line of text 
cat $SCRATCH/pilot_ntm_id.txt | while read fn

do

sample=$fasta/${fn}
echo ${sample}
$salmon_dir/salmon quant 	-i $SCRATCH/gencode_v38_salmon_idx \
							--libType A \
							--skipQuant \
							-o $SCRATCH/salmon_${fn}_quant \
							-1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz \
							-p 4
							
done	
```

## Prepare reference (RSEM + STAR) 

* To optimize, we prepare the reference using **STAR** inside RSEM

```sh

#!/bin/bash
#PBS -l walltime=03:00:00,select=10:ncpus=1:mem=60000mb
#PBS -N RSEM_index_GRCh38.p13
#PBS -A st-bquon-1
#PBS -o output.txt
#PBS -e error.txt

##################################################################
###  RSEM index for alignment
##################################################################

# load perl module
module load perl/5.26.2

# run from scratch/mprieto directory
cd $PBS_O_WORKDIR

# first remove, if it exists, and then create a directory for output
rm -rf index_RSEM
mkdir index_RSEM

# export RSEM software and specify path to STAR and GENCODE
export PATH=/arc/home/mdprieto/software/RSEM-1.3.3/:$PATH
STAR='/home/mdprieto/software/STAR-2.7.8a/source'
GENCODE='/arc/project/st-bquon-1/mprieto/GRCh38.p13/GENCODE'

# specify the reference in .gtf format, argument_1=fasta files, argument_2=prefix $
rsem-prepare-reference  --gtf $GENCODE/gencode.v38.annotation.gtf \
						            --star \
						            --star-path $STAR \
						            -p 8 \					
					              $GENCODE/GRCh38.p13.genome.fa \
						            index_RSEM/GRCh38.p13_RSEM_index \

```

## Quantification step (RSEM + STAR)

* We calculate expression against the previously created reference index
* By calling `STAR` inside `RSEM` the algorithm first aligns and then quantifies, requires a longer walltime but can be done in a single run
* Produces a genome aligned `.bam`

```sh
#!/bin/bash
#PBS -l walltime=90:00:00,select=1:ncpus=10:mem=82000mb
#PBS -N RSEM_counts_NTM_all_STAR
#PBS -o output.txt
#PBS -e error.txt

##################################################################
## RSEM counts
##################################################################

# script to count number of transcripts in our sequencing data
# run from directory /scratch/mprieto and create output directory
cd $PBS_O_WORKDIR
rm -rf counts_NTM
mkdir counts_NTM

# export RSEM & STAR software and load necessary modules (perl, r, gcc)
export PATH=/arc/home/mdprieto/software/RSEM-1.3.3/:$PATH
export PATH=/home/mdprieto/software/STAR-2.7.8a/source:$PATH
module load gcc/5.4.0
module load perl/5.26.2
module load r/3.6.0-py3.7.3

# create PATH variables to different directories
GENCODE='/arc/project/st-bquon-1/mprieto/GRCh38.p13/GENCODE'
STAR='/home/mdprieto/software/STAR-2.7.8a/source'
SCRATCH='/scratch/st-bquon-1/mprieto'
fasta='/scratch/st-bquon-1/mprieto/NTM_fasta'

# txt with sample list must be in scratch directory
cat $SCRATCH/ntm_id.txt | while read fn

do

# sample includes basename of file and path id in the form of CFB2XXX
sample=$fasta/${fn}
id=$(echo ${sample} | sed -r 's/.*(CFB)([0-9]+).*/\12\2/')
echo 'start' ${id}
date

 rsem-calculate-expression -p 10 \
						--paired-end \
						--star \
						--star-gzipped-read-file \
						--star-path $STAR \
						--star-output-genome-bam \
						--strandedness reverse \
						--no-bam-output \
                        ${sample}_R1.fastq.gz  ${sample}_R2.fastq.gz \
                        $SCRATCH/index_RSEM/GRCh38.p13_RSEM_index \
                        $SCRATCH/counts_NTM/${id}
done

# --forward-prob 0 equals ISR libtype / same as --strandedness reverse
# --star uses STAR aligner as part of the counting algorithm
# --star-path specify STAR path too

```

# Alignment QC

## STAR alignment (for picard)

Output a sorted `.bam` to feed directly into picard tools

```sh
#!/bin/bash
#PBS -l walltime=12:00:00,select=1:ncpus=10:mem=60000mb
#PBS -N align_star_ntm
#PBS -m ae
#PBS -o output.txt
#PBS -e error.txt
##################################################################
#### STAR alignment
##################################################################
# run from scratch/mprieto directory
cd /scratch/st-bquon-1/mprieto
# make sure that no such directory exists already
[ -d aligned_STAR ] || mkdir aligned_STAR
# export STAR software
export PATH=/home/mdprieto/software/STAR-2.7.8a/source:$PATH
FASTQ='/scratch/st-bquon-1/mprieto/NTM_fasta'
GENCODE='/arc/project/st-bquon-1/mprieto/GRCh38.p13/GENCODE'
STARindex='/arc/project/st-bquon-1/mprieto/GRCh38.p13/STAR_index'
for sample in $(</scratch/st-bquon-1/mprieto/ntm_id.txt)
do
echo $sample
STAR    --runThreadN 20 \
        --genomeDir $STARindex \
        --sjdbGTFfile $GENCODE/gencode.v38.annotation.gtf \
        --readFilesIn $FASTQ/$sample\_R1.fastq.gz $FASTQ/$sample\_R2.fastq.gz \
        --outFileNamePrefix /scratch/st-bquon-1/mprieto/aligned_STAR/$sample. \
        --sjdbOverhang 100 \
        --readFilesCommand zcat \
        --quantMode TranscriptomeSAM GeneCounts \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM Unsorted SortedByCoordinate
done
```

## Collect RNAseq metrics / fastQC aligned files

Now, we will use `picard_tools` to analyze the quality of the alignment performed by STAR 

* PICARD must be called from java using `java -jar`
* We also run **fastqc** in the aligned `.bam` files in the same script

```sh
#!/bin/bash
#PBS -l walltime=50:00:00,select=1:ncpus=10:mem=60000mb
#PBS -N alignment_QC_picard_fastqc
#PBS -o output.txt
#PBS -e error.txt
##################################################################
#### PICARD CollectRnaSeqMetrics + FastQC of aligned samples
##################################################################

cd /scratch/st-bquon-1/mprieto

# Paths to .bam files, output dir and picard.jar
BAM_DIR='/scratch/st-bquon-1/mprieto/counts_RSEM'
picard='/arc/home/mdprieto/software/picard_tools_2.26'
QC_DIR='/scratch/st-bquon-1/mprieto/QC_ntm/picard_rna_metrics'

# pass all bam files to PICARD command
for bam in `ls $BAM_DIR/*genome.bam`

do 
sample=$(echo ${bam} | egrep -o CFB[0-9]+)
echo ${sample}

java -jar $picard/picard.jar SortSam \
			I=${bam} \
			O=$BAM_DIR/${sample}_sorted.bam \
			SORT_ORDER=coordinate
done

# after we obtain sorted.bam files, we run collectrnaseq metrics
# STRAND defined according to salmon library type results (ISR)
for bam in `ls $BAM_DIR/*sorted.bam`

do
sample=$(echo ${bam} | egrep -o CFB[0-9]+)
echo ${sample}

java -jar $picard/picard.jar CollectRnaSeqMetrics \
			I=${bam} \
			O=$QC_DIR/${sample}_output.RNA_Metrics \
			REF_FLAT=$picard/refFlat.txt \
			STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
			RIBOSOMAL_INTERVALS=$picard/rRNA.interval_list
done

# We then do fastqc on aligned.bam
QC_DIR='/scratch/st-bquon-1/mprieto/QC_ntm'
FASTQC='/arc/home/mdprieto/software/FastQC'
BAM_DIR='/scratch/st-bquon-1/mprieto/counts_RSEM'

# conditional creation of output directory
[ -d QC_ntm/fastqc_aligned ] || mkdir QC_ntm/fastqc_aligned

# Loop genome.bam files to fastqc 
for bam in `ls $BAM_DIR/*genome.bam`

do 
echo running "$bam"
date
echo $FASTQC/fastqc 	-o $QC_DIR/fastqc_aligned \
                -t 8 \
                $bam \
				--noextract
done

```

## Merge all QC using picard

* Change to directory with all QC output
* Run **multiqc** in all subdirectories to produce html report
* Send data to local computer using **sftp** protocol

```sh
export PATH=/arc/home/mdprieto/software/MultiQC:$PATH

cd /scratch/st-bquon-1/mprieto/QC_ntm
multiqc . 
```

* Over-represented sequences continue to be a problem after alignment. Exploratory BLAST shows mostly _hemoglobin subunits_
