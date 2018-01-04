# REMU
REMU is a tool for resolving multi-mappings in read mapping. 
The REMU method uses a Bayesian updating approach for selecting the correct 
alignment for multireads. It can be used as a post-processing step after the
read mapping is done.

## Installation
REMU is developed in pure Python and does not require any installations 
beyond a Python3 distribution. In order to obtain REMU, simply download or clone this
 repository:

 ```commandline
git clone https://github.com/oloomi/remu.git
```

## Getting started
REMU can be used with the following command line instruction:
 
 ```commandline
python3 remu.py -g genome.fna -i input.sam -o output.sam
```
REMU comes with two example files `toy-genome.fna` and `toy-input.sam` for
getting started. The example files can be found in the `example` directory.

## Usage
The input SAM file should contain all mappings. For example, using ```-a``` 
option with BWA and Bowtie2 will report all possible alignments
for a read:

```commandline
bowtie2-build genome.fna genome
bowtie2 -a -x genome -U reads.fq -S all_mappings.sam
```

The multimappings should appear in consecutive lines in the input file,
which is typically the case for most of the alignment softwares. In order to make
sure this requirement is satisfied, the alignments can be sorted based on read name
using the following SAMTools command:
```commandline
samtools sort -n all_mappings.sam -o all_mappings_sorted.sam
```
REMU can then be used on the output of read mapping to resolve the multi-mappings
and select the correct alignment for multi-mapping reads:
```commandline
python3 remu.py -g genome.fna -i all_mappings_sorted.sam -o correct_mappings.sam
```
The output file will contain the alignment for uniquely mapped reads and the
correct alignment for multi-mapping reads and is ready to be used for downstream analysis.

## Options
Please refer to the help section for further options:
```commandline
python3 remu.py -h
```


