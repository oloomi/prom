#!/bin/bash
# Running all scripts

sh run_art_bowtie2_samtools.sh

sh other_methods.sh

cd ..

python3 bayesian-update-main.py

cd bash-scripts

sh run_samtools_corrected_mapping.sh

sh variant_calling.sh

echo "\n=== Everthing done! ===\n"
