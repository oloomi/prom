#!/bin/bash

dir=`pwd`/
script_path="$( cd "$(dirname "$0")" ; pwd -P )/"

# Preparing data
sh ${script_path}get-data.sh

# Running experiments for M. tuberculosis
cd ${dir}mtb/simulated-data/begin-supermax
sh ${script_path}pipeline.sh -s -b 150

cd ${dir}mtb/simulated-data/middle-supermax
sh ${script_path}pipeline.sh -s -m 150

cd ${dir}mtb/real-data/back-mutate
sh ${script_path}pipeline.sh -r -b 150

# Running experiments for E. coli
cd ${dir}ecoli/simulated-data/begin-supermax
sh ${script_path}pipeline.sh -s -b 150

cd ${dir}ecoli/simulated-data/middle-supermax
sh ${script_path}pipeline.sh -s -m 150

cd ${dir}ecoli/real-data/back-mutate
sh ${script_path}pipeline.sh -r -b 100

# Running experiments for S. cerevisiae
cd ${dir}yeast/simulated-data/begin-supermax
sh ${script_path}pipeline.sh -s -b 150

cd ${dir}yeast/simulated-data/middle-supermax
sh ${script_path}pipeline.sh -s -m 150

cd ${dir}yeast/real-data/back-mutate
sh ${script_path}pipeline.sh -r -b 150