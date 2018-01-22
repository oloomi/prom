#!/bin/bash
# Running all scripts
# Run this script in the specific experiment folder

# sh all.sh -s -b
# sh all.sh -s -m
# sh all.sh -r

dir="/mnt/remu/evaluation/"

# Simulated data
if [ "$1" = "-s" ]; then
  # Mutations in the beginning of repeats
  if [ "$2" = "-b" ]; then
    python3 $dir"evaluation.py" -cb
  # Mutations in the middle of repeats
  elif [ "$2" = "-m" ]; then
    python3 $dir"evaluation.py" -cm
  else
    echo "Type of simulation required!"
  fi
# Real data
elif [ "$1" = "-r" ]; then
  python3 $dir"evaluation.py" -b
else
  echo "Invalid argument in all.sh!"
fi

sh $dir"bash-scripts/read-mapping.sh" $1
sh $dir"bash-scripts/multimapping-resolution.sh" $1
sh $dir"bash-scripts/variant-calling.sh" $1
python3 $dir"evaluation.py" -e

echo "\n=== Everthing done! ===\n"


