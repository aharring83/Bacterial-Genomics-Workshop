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
rm -rf test
```

### Quality Control Assessment and Trimming of Sequencing Reads

A FASTQ file is a text file that contains the sequence data.
[![](https://drive5.com/usearch/manual/fastq_fig.jpg)

The fastq files are in a folder "Data", located at:
These are the raw fastq files that have not been trimmed of low quality data and sequencing adapters/indexes.
We are going to trim and assess the quality of the trimmed reads.


