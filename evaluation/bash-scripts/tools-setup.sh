#!/usr/bin/env bash

sudo apt-get update
sudo apt-get install gcc
sudo apt install make
sudo apt install vim
sudo apt install unzip
sudo apt install git
sudo apt install g++

# SAMtools dependencies
sudo apt-get install libncurses5-dev
sudo apt-get install zlib1g-dev
sudo apt-get install libbz2-dev

# SAMtools
sudo apt-get install liblzma-dev    mkdir /mnt/samtools
sudo apt-get install curl           wget "https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2"
tar xvjf samtools-1.6.tar.bz2
cd samtools-1.6
./configure --prefix=/mnt/tools/samtools
make
make install

# BCFtools
wget "https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2"

# Bowtie2
wget "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip"
unzip bowtie2-2.3.3.1-linux-x86_64.zip
mv bowtie2-2.3.3.1-linux-x86_64 bowtie2

# BWA
wget "https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2"
tar xjvf
make

# Freebayes
git clone --recursive https://github.com/ekg/freebayes.git
make
sudo make install

# Path variable
export PATH=$PATH:/mnt/tools/samtools/bin
export PATH=$PATH:/mnt/tools/bcftools/bin
export PATH=$PATH:/mnt/tools/bowtie2

source ~/.profile

