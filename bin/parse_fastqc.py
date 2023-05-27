#!/usr/bin/env python

import re
import subprocess

import sys

'''
Script for extracting modules from fastqc.txt output and doing visualizations
'''


def format(string: str) -> str:
    string = string.replace('>>', '').replace(' ', '_')
    string = re.sub('\t.*', '', string)
    return string.strip() + '.txt'


def get_num(string: str) -> int:
    return int(re.sub('[^0-9]+', '', string))


directory = sys.argv[1]
dir_name: str = directory[:directory.find('_')]
locs = subprocess.run([f'cat {directory}/fastqc_data.txt | grep ">>" -n'],
                      shell=True,
                      stdout=subprocess.PIPE)
stdout = str(locs.stdout).split('\\n')[:-1]
modules = [get_num(num) for num in stdout if '_' not in num]
ends = [get_num(num) for num in stdout if 'END' in num]
with open(f'{directory}/fastqc_data.txt') as f:
    to_parse = f.readlines()
    for strt, end in zip(modules, ends):
        module = format(to_parse[strt-1])
        formatted = open(f'{dir_name}_{module}', 'w')
        formatted.write(''.join(to_parse[strt:end-1]))
