#!/bin/bash

dir=`pwd`/
script_path="$( cd "$(dirname "$0")" ; pwd -P )/"

# Running experiments
cd ${dir}mtb/simulated-data/begin-supermax
sh ${script_path}bwa-mmr-fix-pipeline.sh -s 150

cd ${dir}mtb/simulated-data/middle-supermax
sh ${script_path}bwa-mmr-fix-pipeline.sh -s 150

cd ${dir}mtb/real-data/back-mutate
sh ${script_path}bwa-mmr-fix-pipeline.sh -r 150

# Running experiments
cd ${dir}ecoli/simulated-data/begin-supermax
sh ${script_path}bwa-mmr-fix-pipeline.sh -s 150

cd ${dir}ecoli/simulated-data/middle-supermax
sh ${script_path}bwa-mmr-fix-pipeline.sh -s 150

cd ${dir}ecoli/real-data/back-mutate
sh ${script_path}bwa-mmr-fix-pipeline.sh -r 150

# Running experiments
cd ${dir}yeast/simulated-data/begin-supermax
sh ${script_path}bwa-mmr-fix-pipeline.sh -s 150

cd ${dir}yeast/simulated-data/middle-supermax
sh ${script_path}bwa-mmr-fix-pipeline.sh -s 150

cd ${dir}yeast/real-data/back-mutate
sh ${script_path}bwa-mmr-fix-pipeline.sh -r 100