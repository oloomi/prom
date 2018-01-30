#!/bin/bash
# Running all scripts
# Run this script in the specific experiment folder

# sh experiment.sh -s -b 150
# sh experiment.sh -s -m 150
# sh experiment.sh -r -b 100

script_path="$( cd "$(dirname "$0")" ; pwd -P )/"

# Simulated data
if [ "$1" = "-s" ]; then
  # Mutations in the beginning of repeats
  if [ "$2" = "-b" ]; then
    python3 ${script_path}"evaluation.py" -sb $3
  # Mutations in the middle of repeats
  elif [ "$2" = "-m" ]; then
    python3 ${script_path}"evaluation.py" -sm $3
  else
    echo "Type of simulation required!"
  fi
# Real data
elif [ "$1" = "-r" ]; then
  if [ "$2" = "-b" ]; then
    python3 ${script_path}"evaluation.py" -rb $3
  fi
else
  echo "Invalid argument in all.sh!"
fi

sh ${script_path}"read-mapping.sh" $1 $3
sh ${script_path}"multimapping-resolution.sh" $1 $3
sh ${script_path}"variant-calling.sh" $1
python3 ${script_path}"evaluation.py" -e

echo "\n=== Experiment completed! ===\n"


