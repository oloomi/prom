
"""
Copyright (C) 2017 S Mohammad H Oloomi (smh.oloomi@gmail.com)

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""
import bayesian_update as bu
import argparse
import textwrap
import os
import psutil

parser = argparse.ArgumentParser(usage='\n\tpython3 remu.py -i <input.sam> -r <reference.fa> [-o <output.sam>]\n' +
                                       'help:\n\tpython3 remu.py -h',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                                 REMU: REsolving MUltimappings in Read Mapping
                                 by Mohammad Oloomi (smh.oloomi@gmail.com)\n'''),
                                 epilog='')
parser.add_argument('-i', '--input', required=True, help='Input SAM file that contains all read mappings')
parser.add_argument('-r', '--reference', required=True, help='Reference genome in FASTA format')
parser.add_argument('-o', '--output', default='output.sam',
                    help='Output SAM file that contains correct mappings (default: output.sam)')
# Software version
parser.add_argument('-v', '--version', action='version', version='REMU 1.0.0', help='Software version')
args = parser.parse_args()

if args.input and args.reference:
    print(parser.description)
    bu.bayesian_update(args.reference, args.input, args.output)

    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    print("Memory:", mem, "MB")


# -r "./data/genomes/toy-genome.fna" -i "./read-mapping/toy-genome-mutated/toy-wg-mutated-se-mapping-report-all.sam" -o "./read-mapping/toy-genome-mutated/test.sam"

# python3 remu.py -r "../data/NC_000962.fna" -i "../data/toy-report-all.sam" -o "../data/remu-toy-report-all.sam"