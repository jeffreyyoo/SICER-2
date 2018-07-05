# SICER 2.0

## Introduction
*to be filled out*

For more information about the original SICER algorithm, please see

>> “*A clustering approach for identification of enriched domains from histone modification
>> ChIP-Seq data” Chongzhi Zang, Dustin E. Schones, Chen Zeng, Kairong Cui, Keji Zhao, and
>> Weiqun Peng*, **Bioinformatics** 25, 1952 - 1958 (2009)

In addition, we present an alternative algorithm for *insert description of recognicer*

For more information about the RECOGNICER algoirthm, please see

>> *insert citation for RECOGNICER paper*

## Installation
### Requirements
#### Python Version
Unlike the original version of SICER, SICER 2.0 runs in **Python 3**.
Please use Python 3 to install and run SICER 2.0. The program should work with most of versions of python3.

#### Libraries
Numpy and Scipy are required to run SICER. Please have these installed before installing SICER.
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

Note that SICER 2.0 has two algorithms for the users can use to find candidats for significant islands. Sub-commands `SICER` and `RECOGNICER` are used to decide which one to use.

### SICER
To use the SICER algorithm, type the subcommand `SICER` after the `sicer` command (i.e. type `sicer SICER` in terminal).

#### Arguments
##### -t/--treatment_file (Required)
The file must either be in BED or BAM format (note that you need *bedtools* installed to directly enter BAM files).
The file name can either just be the name or the full path of the file. 

##### -c/--control_file (Optional)
Like the treatment file, control file must be in BED or BAM format. However, control library input is optional.

##### -s/--species (Required)
ex) `-s hg38`

##### -rt/--redundancy_threshold (Optional)
The number of copies of indentical reads allowed in a library. Default value is 1

##### -w/--window_size (Optional)
Resolution of SICER. Default value is 200 (bp)

##### -f/--fragment_size (Optional)
The amount of shift from the beginning of a read to the center of the DNA fragment represented by the read. 
Default value is 150 (bp).

##### -egf/--effective_genome_faction (Optional)
Effective genome as fraction of the genome size. Default value is 0.74.

##### -fdr/--false_discovery_rate (Optional)
Remove all islands with an false-discovery-rate below cutoff. Default value is 0.01.

##### -g/--gap_size (Optional)
The minimum length of a "gap" such that neighboring window is an "island." 
Please note that this value must be a multiple of the window size.
Default value is 600 (bp).

##### -e/--e_value (Optional)
E-value. Requires user input when no control library is provided. Default value is 1000

##### -i/--input_directory (Optional)
Path of the directory in which input files are located in. Not required when running SICER in the same directory as input files.
Default is the current working directory.

##### -o/--output_directory (Optional)
Path of the directory in which results will be stored. Default output directory is the current working directory.

##### -opt_o/--optional_output (Optional)
Additional Outputs: Enter "True" or "1" to have SICER produce a BED file of treatment reads filtered by significant islands and WIG file of filtered reads binned into windows.

### RECOGNICER

To use the RECOGNICER algorithm type the subcommand `RECOGNICER` (e.g. `sicer RECOGNICER`)

#### Arguments
All of the arguments for RECOGNICER are identical to those of SICER except for `gap_size` and `e_value`. 
Instead of these two arguments, RECOGNICER has two arguments called `step_size` and `step_score`.

##### -s_size/--step_size (Optional)
The number of windows in one graining unit. Default value is 3.

##### -s_score/--step_score (Optional)
The minimum number of positive elements in the graining unit to call the unit positive. Default value is 2.

### Differential Peak Calling
The command for differential peak calling is `sicer_df`. So for example, if you wish to run differential peak calling
using the RECOGNICER algorithm, you would type `sicer_df RECOGNICER`.

#### Arguments
Most of the arguments for both SICER and RECOGNICER differential peak calling are identical to those of the regular peak callings except for the following arguments specified below. 

Also, differential peak calling has one additioanl argument called `----false_discovery_rate_df`

##### -t/--treatment_file (Required)
Two files must be given as input. The first file must be the knockout (KO) file and the second file must be the wild-type (WT) file.
Both files must either be in BED or BAM format. 

##### -c/--control_file (Optional)
While optional, two files must be given as input if you decide to provide the input. The first file must be the control library corresponding to the knockout (KO) treatment file and the second file must be the control library corresponding to the wild-type (WT) treatment file. Both files must either be in BED or BAM format.

##### -fdr_df/--false_discovery_rate_df (Optional)
Cutoff for identification of significant changes been wild-type library and knockout library. Default value is 0.01. 


## Example Calls
1. Calling SICER with a control library. 
*Default parameters are explicitly entered for the sake of demonstration.*

`sicer SICER -t treatment.bed -c control.bed -s hg38 -w 200 -rt 1 -f 150 -egf 0.74 -fdr 0.01 -g 600 -e 1000`

2. Calling SICER without a control library

`sicer SICER -t treatment.bed -s hg38`

3. Calling SICER with control libraries for differential peak calling.

`sicer_df SICER -t treatment1.bed treatment2.bed -c control1.bed control2.bed -s hg38`

4. Calling SICER without control libraries for differential peak calling.

`sicer_df SICER -t treatment1.bed treatment2.bed -s hg38`

Replace the second word "SICER" with "RECOGNICER" to use RECOGNICER algorithm.

## Questions?
For technical questions or issues, feel free to contact Jin Yong (Jeffrey) Yoo at jy2ma@virginia.edu.
For questions about the methodology, please contact Chongzhi Zang, PhD. at zang@virginia.edu. 
