#!/usr/bin/python3

"""
Copyright (C) 2017 S Mohammad H Oloomi (smh.oloomi@gmail.com)

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""
import bayesian_update as bu
import argparse
import textwrap


def float_prob(num_str):
    num = float(num_str)
    if not (0 <= num <= 1):
        raise argparse.ArgumentTypeError("{} is not in range [0, 1]".format(num))
    return num

parser = argparse.ArgumentParser(usage='\n\tpython3 remu.py -i <input.sam> -r <reference.fa> [-o <output.sam>]\n' +
                                       'help:\n\tpython3 remu.py -h',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                                 REMU: REsolving MUltimappings in Read Mapping
                                 by Mohammad Oloomi (smh.oloomi@gmail.com)\n'''),
                                 epilog='')
parser.add_argument('-i', '--input', required=True, help='Input SAM file that contains all read mappings')
parser.add_argument('-g', '--genome', required=True, help='Reference genome in FASTA format')
parser.add_argument('-o', '--output', default='output.sam',
                    help='Output SAM file that contains correct mappings (default: output.sam)')
parser.add_argument('-r', '--runs', default=10, type=int, help='Number of runs (default: 10)')
parser.add_argument('-t', '--threshold', default=0.02, type=float_prob,
                    help='Threshold for probability difference between equally good mappings (default: 0.02)')
# Software version
parser.add_argument('-v', '--version', action='version', version='REMU 1.0.0', help='Software version')
args = parser.parse_args()

if args.input and args.genome:
    print(parser.description)
    bu.bayesian_update(args.genome, args.input, args.output, args.runs, args.threshold)
