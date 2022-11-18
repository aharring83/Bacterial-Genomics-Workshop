# Bacterial-Genomics-Workshop
- Welcome to the bacterial genomics workshop where we will learn how to process sequencing data, assemble genomes and build phylogenetic trees.
- We will being using a remote server to do the analysis but the procedures can be performed on a personal laptop with conda installed. 
(- Instructions can on installing conda can be found here: https://github.com/rltillett/conda_notes


------------
# Requirements:
Laptop (Mac or PC)
If you have a PC laptop, please install gitbash or putty prior to the workshop.

------------

### Command Line Interface (CLI)
Most bioinformatic tools use the CLI. The CLI is a way of interacting with a computer's operating system using text representation of files, folders, etc. This differs from the commonly used graphical user interface (GUI), where we can see images linked to files, folders, etc. In order to access the CLI, we have to invoke a program, shell, which is a CLI program that allows us to control the computer by typing instructions with a keyboard.

On a Mac or Linux machine, you can access a shell through a program called “Terminal”, which is already available on your computer.
On a PC/Windows, you’ll need to download powershell or gitbash to access the shell.

To save time, we will be working on a remote server, provided by CyVerse (located in Arizona), where all the necessary data and software is available. 

Here is a cheatsheet to help you with using text commands to interact with the CLI.

[![](https://rumorscity.com/wp-content/uploads/2014/08/10-Linux-Unix-Command-Cheat-Sheet-021.jpg)](http://https://rumorscity.com/wp-content/uploads/2014/08/10-Linux-Unix-Command-Cheat-Sheet-021.jpg)

Lets begin using the CLI!!!!

------------ 

First we need to connect to the remote server using secure shell (ssh). For Macs, open a terminal. For PC, open a terminal using gitbash or putty. Type the following and press enter:
```
ssh {username}@10.238.4.152
```
Press enter and you will be prompted to enter your password. Enter your password, press enter. You will be asked a question about connecting to the server for the first time. Type "Yes" and press Enter.

Congratulations! You have now used the CLI to connect to a remote server.
Now we will test a few commands to get more familar with the CLI. Try the following commands:
```
ls
```
```
pwd
```
```
mkdir test
```
```
cd test
```
```
pwd
```
```
cd ..
```
```
ls 
```
```
rm -rf test
```
```
ls
```

### Activating conda and creating an environment
Exectue the following commands:
```
tmux
```
```
source ../miniconda/bin/activate
```
Lets activate our environment for the workshop which has all the packages we need for the workshop installed already 
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
```
conda activate workshop
```

Now we are readying to do some bioinformatics!!!

### Quality Control Assessment and Trimming of Sequencing Reads

A FASTQ file is a text file that contains the sequence data.

[![](https://external-content.duckduckgo.com/iu/?u=https%3A%2F%2Fi1.wp.com%2Fgencoded.com%2Fwp-content%2Fuploads%2F2020%2F05%2Ffastq_format_explained-3.png%3Fssl%3D1&f=1&nofb=1&ipt=48754aaa7290bcaf4d9aa3c070fde25c449e2451a8d775d00849567eb23ba5b9&ipo=images)

The fastq files are in a folder "data", located one directory above your home directory. Lets copy it to your directory.
```
cp -r ../data .
```
These are the raw fastq files that have not been trimmed of low quality data and sequencing adapters/indexes.
We are going to trim and assess the quality of the trimmed reads.
We will be using fastp to trim the reads, FastQC and MultiQC to assess the reads quality.
Since we have a few samples to process, it would be easier if we wrote a bash script to execute the commands.
We will be using a text editor, nano, to write this script, but any text editor will work.
Enter the commands:
```
nano
```
Paste this code:
```
#!/bin/bash
mkdir -p trim
mkdir -p QC
for i in *_R1.fastq.gz
do
j=${i/_R1.fastq.gz/_R2.fastq.gz}
fastp -i ${i} -I ${j} -o trim/${i/_R1.fastq.gz/_trim_R1.fastq.gz} -O trim/${j/_R2.fastq.gz/_trim_R2.fastq.gz} -w 16
wait
fastqc -t 32 -o QC trim/*.gz
done
```
Press Ctrl + X. Type Yes, name the file "trimQC.sh"
```
ls
```
```
chmod u+x trimQC.sh
```
```
ls
```
We have create a shell script that is executable!
Lets execute the script.
```
./trimQC.sh
 ```
 
 After the script is running, we should have 2 new folders (trim and QC) created. 
 Lets view the quality of our trimmed reads. We need to combine all the results into one file using MultiQC. Run this command:
 ```
 multiqc QC
 ```
 Since our remote server is headless, we cannot view the files on the server, so we need to transfer the files to our personal computer/laptop to view.
 Open a new terminal and enter:
 ```
 scp {username}@10.238.4.152:/home/{username}/data/multiqc_report.html .
 ```
 Open the downloaded html file in a web browser.
 
 ### Reference-based genome assembly
[![](https://i.ibb.co/jLMwQt0/Picture1.png)

We are going to use Burrows-Wheeler aligner to align our reads and ivar to generate a consensus genome. We will then use bcftools to generate a variant calling file (vcf), which contains single nucleotide polymorphisms of each sample compared to the reference genome. Finally, we will convert the vcf files to fasta files, perform a multiple sequence alignment and end the workshop with phylogenetic analysis.

First, we need to create an index of our reference genome.
Run the following commands in the ref folder, which is inside the data folder:
```
bwa index ref.fasta ref.fasta
```
```
ls
```
Now we will align the reads to the indexed reference genome. Make sure you are in the data directory and not the ref directory!
```
nano
```
Paste this code:
```
#!/bin/bash
#$1 folder containing trimmed reads
for i in $1/*_trim_R1.fastq.gz
do
prefix=$(basename ${i/_trim_R1.fastq.gz})
j=${i/_trim_R1.fastq.gz/_trim_R2.fastq.gz}
mkdir -p sam
mkdir -p bam
bwa mem -t 16 -o sam/${prefix}.sam /home/ref/ref.fasta ${i} ${j}
samtools view -h -b sam/${prefix}.sam| samtools sort -@16 -o bam/${prefix}.bam
samtools mpileup -aa -A -d 0 -Q 0 bam/${prefix}.bam} | ivar consensus -p ${prefix}
lofreq call -f /home/ref/ref.fasta -o ${prefix}.vcf bam/${prefix}.bam
wait
bedtools getfasta -fi /home/ref/ref.fasta -bed ${prefix}.vcf -fo ${prefix}_vcf.fasta
grep -v '>' ${prefix}_vcf.fasta | tr -d  '\n' > ${prefix}_vcf.fasta}
sed -e '1i\>' ${prefix}_vcf.fasta > ${prefix}_snp.fasta
sed -i "s/^>.*/&${prefix}/" ${prefix}_snp.fasta

done



#!/bin/bash
# folder containing trimmed reads is $1
# output folder is $2
mkdir -p $2/sam
mkdir -p $2/bam
mkdir -p $2/vcf
mkdir -p $2/fasta
for i in $1/*_trim_R1.fastq.gz
do
prefix=$(basename ${i/_trim_R1.fastq.gz})
j=${i/_trim_R1.fastq.gz/_trim_R2.fastq.gz}
sam=$2/sam/${prefix}.sam
bam=$2/bam/${prefix}.bam
vcf=$2/vcf/${prefix}.vcf
fasta=$2/fasta/${prefix}_vcf.fasta
fastatree=$2/fasta/${prefix}_tree.fasta
fastasnp=$2/fasta/${prefix}_snp.fasta
bwa mem -t 16 -o ${sam} /home/ref/ref.fasta ${i} ${j}
samtools view -h -b ${sam}| samtools sort -@16 -o ${bam}
#samtools mpileup -aa -A -d 0 -Q 0 ${bam} | ivar consensus -p ${prefix}
lofreq call -f /home/ref/ref.fasta -o ${vcf} ${bam}
wait
bedtools getfasta -fi /home/ref/ref.fasta -bed ${vcf} -fo ${fasta}
grep -v '>' ${fasta} | tr -d  '\n' > ${fastatree}
sed -e '1i\>' ${fastatree} > ${fastasnp}
sed -i "s/^>.*/&${prefix}/" ${fastasnp}
done
```
Press Ctrl + X. Type Yes, name the file "analysis.sh"
```
ls
```
```
chmod u+x analysis.sh
```
```
ls
```
```
./analysis.sh trim
 ```

 ### Multiple Sequence Alignment and Phylogenetic Analysis based on SNPs.
 
We are going to use MAFFT to do the multiple sequence alignment and fasttree to generate a phylogenetic tree.
Run the following commands:
```
cat *_snp.fasta > all.fasta
```
```
mafft --thread 44 all.fasta > align.fasta
```
```
FastTree -nt align.fasta > align.tree
```
```
Open a new terminal and transfer the align.tree file to your personal computer
Go to http://etetoolkit.org/treeview/
Upload your align.tree file and press view tree.
Congratulations, we have generated a phylogenetic tree based on the individual Legionella species SNVs.
Scientist use these trees to trace origins and spread of outbreaks.
From this workshop, we can use the tree to see which Legionella strains our samples are related to.
Thank you for participating in the workshop.
```




