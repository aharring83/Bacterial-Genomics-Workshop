# Bacterial-Genomics-Workshop
- Welcome to the bacterial genomics workshop where we will learn how to process sequencing data, assemble genomes and build phylogenetic trees.
- We will being using CyVerse to conduct this workshop but the procedures can be performed on a personal laptop with conda installed. 
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

First we need to connect to the CyVerse server which offers cloud-computing. We need to launch an instance of a pre-configured virtual machine. CyVerse uses a program, Atmosphere, to launch instances. We will be starting our atmospheres (instances) and connecting to the remote server by ssh (secure shell).
Login into CyVerse, then launch Atmosphere and login.  
