#!/usr/bin/env python3
import sys
import subprocess

import argparse

parser = argparse.ArgumentParser(description='Description.')
parser.add_argument('--in', dest='input', required=True, help='Flagstat input file')
parser.add_argument('--out', dest='output', required=True, help='Percentage output file.')
args = parser.parse_args()


def read_mapping_info(flagstats_file):
    flagstat_info = {}
    with open(flagstats_file, 'r') as reader:
        for line in reader:
            if 'total' in line:
                flagstat_info['total'] = int(line.split(' ')[0])
            if 'mapped (' in line:
                flagstat_info['mapped'] = int(line.split(' ')[0])
    return flagstat_info

def calculate_percentage_mapped(flagstat_info):
    mapped = flagstat_info['mapped']
    total = flagstat_info['total']
    if total == 0:
        return 0.0
    else:
        return round( mapped / total * 100, 2 )


infos_from_flagstat = read_mapping_info(args.input)
percentage = str(calculate_percentage_mapped(infos_from_flagstat))

with open(args.output, 'w') as writer:
    writer.write(percentage)