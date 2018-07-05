# SICER 2.0

## Introduction
??

For more information about the original SICER algorithm, please see

>> “*A clustering approach for identification of enriched domains from histone modification
>> ChIP-Seq data” Chongzhi Zang, Dustin E. Schones, Chen Zeng, Kairong Cui, Keji Zhao, and
>> Weiqun Peng*, **Bioinformatics** 25, 1952 - 1958 (2009)

In addition, we present a alternative algorithm for *insert description of recognicer*

For more information about the RECOGNICER algoirthm, please see

>> *insert citation for RECOGNICER paper*

## Installation
### Requirements
#### Python Version
Unlike the original version of SICER, SICER 2.0 runs in **Python 3**.
Please use Python 3 to install and run SICER 2.0. The program should work with most of versions of python3.

#### Libraries
Numpy and Scipyare required to run SICER. Please have these installed before installing SICER.
This can be done by simply typing `pip install numpy scipy` under command line (if python2.7 is your default python version, use `pip3`).

#### C Compiler
C compiler is required to compile C codes that are part of the SICER package. This also means that python header files (e.g. Python.h) are needed. For Linux users, make sure to have python-dev installed. For Mac OS X users, it is recommneded that you install Xcode.

#### BedTools
Lastly, if you would like to directly pass BAM files as input files for SICER, you need to have *bedtools* installed. Please refer to this [link](http://bedtools.readthedocs.io/en/latest/) for more details on installing bedtools. This is not required if you will intend to only pass BED files as input files.

### Installation
To install SICER through PyPI, simply open the terminal and type `pip install SICER` (or pip3 if python2.7 is your default python). 
To update SICER, you can type in `pip install -U SICER`

## Using SICER
The termianl command to run SICER is `sicer`. 
### Choosing between SICER and RECOGNICER
SICER 2.0 has two algorithms for the users can use to find candidats for significant islands.

To use SICER algorithm, type the subcommand `SICER` after the `sicer` command (i.e. type `sicer SICER` in terminal)

To use RECOGNICER algorithm type the subcommand `RECOGNICER` (e.g. `sicer RECOGNICER`)
There are multiple modes when it comes to running SICER. 

### 1. Running SICER with control library


### 2. Running SICER to Identify Differentially Enriched Regions

