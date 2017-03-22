#!/bin/bash

#art="E:\BioinfTools\art\art_illumina"
export PATH=$PATH:/home/mohammad/Applications/art_bin_MountRainier/
art="art_illumina"

mtb_genome="./mtb-genome-extract-mutated.fna"

#Simulation of single-end reads of 150 bp with coverage 10; with max number of indels = 0 (-k 0)
$art -ss HS25 -sam -i $mtb_genome -l 150 -f 10 -k 0 -rs 12345 -o mtb-single-end

