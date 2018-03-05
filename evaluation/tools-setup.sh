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
sudo apt-get install liblzma-dev
sudo apt-get install curl
mkdir /mnt/samtools
wget "https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2"
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

# ART
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz
tar -xzvf artbinmountrainier20160605linux64tgz.tgz

# Vmatch
wget http://www.vmatch.de/distributions/vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
tar -xzvf vmatch-2.3.0-Linux_x86_64-64bit.tar.gz

# MMR
git clone https://github.com/ratschlab/mmr.git
make

# REMU
git clone https://github.com/oloomi/remu.git
chmod +x remu.py

# SRA toolkit fastq-dump
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz
mv sratoolkit.2.8.2-1-ubuntu64 sratoolkit

# Fastx toolkit
sudo apt install fastx-toolkit

# seqtk
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make

# Path variable
export PATH=$PATH:/home/ubuntu/tools/samtools/bin
export PATH=$PATH:/home/ubuntu/tools/bcftools/bin
export PATH=$PATH:/home/ubuntu/tools/bowtie2
export PATH=$PATH:/home/ubuntu/tools/art_bin_MountRainier
export PATH=$PATH:/home/ubuntu/tools/vmatch-2.3.0-Linux_x86_64-64bit
export PATH=$PATH:/home/ubuntu/tools/mmr
export PATH=$PATH:/home/ubuntu/tools/sratoolkit/bin
source ~/.profile

