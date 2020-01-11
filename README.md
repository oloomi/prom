# PROM
PROM (Probabilistic Resolution Of Multi-mappings) is a tool for resolving multi-mappings in read mapping. 
The PROM method uses a probabilistic approach inspired by Bayesian updating for selecting the correct 
alignment for multi-reads. It can be used as a post-processing step after the
read mapping is done.

## Installation
PROM is developed in pure Python and does not require any installations 
beyond a Python3 distribution. In order to obtain PROM, simply download or clone this
 repository:

 ```commandline
git clone https://github.com/oloomi/prom.git
```

## Getting started
PROM can be used with the following command line instruction:
 
 ```commandline
python3 prom.py -g genome.fna -i input.sam -o output.sam
```
PROM comes with two example files `toy-genome.fna` and `toy-input.sam` for
getting started. The example files can be found in the `example` directory.

## Usage
The input SAM file should contain all mappings. For example, using ```-a``` 
option with BWA and Bowtie2 will report all possible alignments
for a read:

```commandline
bowtie2-build genome.fna genome
bowtie2 -a -x genome -U reads.fq -S all_mappings.sam
```

The multi-mappings should appear in consecutive lines in the input file,
which is typically the case for most of the alignment softwares. In order to make
sure this requirement is satisfied, the alignments can be sorted based on read name
using the following SAMTools command:
```commandline
samtools sort -n all_mappings.sam -o all_mappings_sorted.sam
```
PROM can then be used on the output of read mapping to resolve the multi-mappings
and select the correct alignment for multi-mapping reads:
```commandline
python3 prom.py -g genome.fna -i all_mappings_sorted.sam -o correct_mappings.sam
```
The output file will contain the alignment for uniquely mapped reads and the
correct alignment for multi-mapping reads and is ready to be used for downstream analysis.

## Options
Please refer to the help section for further options:
```commandline
python3 prom.py -h
```

## Citation
Oloomi, S. M. H. (2018). *The impact of multi-mappings in short read mapping* (Doctoral dissertation).

@phdthesis{oloomi2018impact,
  title={The impact of multi-mappings in short read mapping},
  author={Oloomi, Seyed Mohammad Hossein},
  year={2018}
} 

Please feel free to contact me (Mohammad) for any enquiries: [smh.oloomi@gmail.com](mailto:smh.oloomi@gmail.com)
 
 


