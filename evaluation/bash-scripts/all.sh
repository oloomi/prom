#!/bin/bash
# Running all scripts

python3 /mnt/remu/evaluation/evaluation.py -c
sh /mnt/remu/evaluation/bash-scripts/read-mapping.sh
sh /mnt/remu/evaluation/bash-scripts/multimapping-resolution.sh
sh /mnt/remu/evaluation/bash-scripts/variant-calling.sh
python3 /mnt/remu/evaluation/evaluation.py -e

echo "\n=== Everthing done! ===\n"

#├── genome-ref
#│   └── repeats
#├── real-data
#└── simulated-data
#    └── supermax-100-140
#        ├── genome-mutated
#        ├── mappings
#        │   ├── bowtie
#        │   └── bwa
#        ├── reads
#        ├── results
#        └── variants
