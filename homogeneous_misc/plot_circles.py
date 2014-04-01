from __future__ import print_function
import argparse
import os
import sys
import plotdata
import math
from matplotlib import pyplot

"""Read output of runs and summarise per group"""

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', default='filename')
opt = parser.parse_args()

clone = []

lines = [line.strip() for line in open(opt.filename)]
for l in lines:
    c_line = l.split(' ')
    c = [float(c_line[0]),float(c_line[1]),float(c_line[2])]
    clone.append(c)
    

plotdata.mutation_v_proliferation(clone,opt.filename+"plot","circles",1)
