from __future__ import print_function
from matplotlib.patches import Polygon
from prettyplotlib import plt
import prettyplotlib as ppl
import numpy as np
import os
import math
import argparse
from prettyplotlib import mpl


parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename',default='filename')
parser.add_argument('-o','--output',default='output')
opt = parser.parse_args()

filename = opt.filename

data = [[],[],[],[],[],[],[],[],[]]

lines = [line.strip() for line in open(filename)]
locat = 0
for l in lines:
    if l == '':
        print("ignore empty")
    else:
        if l == 'END':
            locat+=1
            print("increasing")
        else:
            n = l.split(' ')
            #print(n)
            value = float(n[0])  #proliferation
            value_mut = float(n[1])
            count = int(n[2])
            group = n[2]
            if value < 1:
                for x in range(count):
                    data[locat].append(value)
        

print("data")
for i in data:
    print(i[:10])
print("data")

#np.random.seed(10)

#data = np.random.randn(8, 4)

u_labels = ['1e-7','1e-6','1e-5','1e-4','0.001','0.005','0.01','0.05','0.1']

fig, ax = plt.subplots()

ax.set_xticklabels(u_labels)
ax.set_xticks(u_labels)
ax.set_xlabel('Initial Mutation Rates',fontsize=9)
ax.set_ylabel('End Proliferation Rates',fontsize=9)
i = 0
print(ax.get_xticklabels())
ppl.boxplot(ax, data) #,xticklabels=u_labels)
plt.title("Distribution of End Proliferation Rates by Initial Mutation Rate",fontsize=9)
fig.savefig(opt.output+'_boxplot.png')
print(opt.output+"DONE")
#ppl.hist(ax,data)
#fig.savefig('histogram_prettyplotlib_default.png')
