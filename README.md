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
Unlike the original version of SICER, SICER 2.0 runs in python3.

Please use **python3** to install and run SICER 2.0. The program should work with most of versions of python3.

Also, numpy and scipy are required to run SICER. Please have these installed before installing SICER. 
This can be done by simply typing in `pip install numpy` and `pip install scipy` (note that numpy must be installed before installing scipy). If python2.7 is your default python version, remember to call `pip3` to install these.

GCC is required to compile .c codes in SICER 2.0 package, and python header files are needed. If you are using Mac OSX, I recommend you install Xcode; if you are using Linux, you need to make sure python-dev is installed.

### Easy Installation Through PyPI
To install SICER through PyPI, simply open the terminal and type in `pip install SICER` (or pip3 if python2.7 is your default python). 

To update SICER, you can type in `pip install -U SICER`



