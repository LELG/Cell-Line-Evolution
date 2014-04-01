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
parser.add_argument('-m','--filemid',default='filename')
parser.add_argument('-f','--filename',default='filename')
parser.add_argument('-o','--output',default='output')
parser.add_argument('-t','--ty',default='output')
opt = parser.parse_args()

filename = opt.filename

data = [[],[],[],[],[],[],[],[],[]]
midd = [[],[],[],[],[],[],[],[],[]]


lines = [line.strip() for line in open(filename)]
locat = 0
for l in lines:
    if l == '':
        print('skip header')
    else:
        #if l == '' or l == 'END':
        if l == 'END':
            locat+=1
            print("next set")
        else:
            n = l.split(' ')
            #value = float(n[j])  
            value = float(n[0])  #mutation rate
            if value < 1:
                data[locat].append(value)
    

for i in data:
    print(i[:10])

print("HALF WY THWRE")

lines = [line.strip() for line in open(opt.filemid)]
locat = 0
for l in lines:
    if l == '':
        print('skip header')
    else:
        #if l == '' or l == 'END':
        if l == 'END':
            locat+=1
            print("next set")
        else:
            n = l.split(' ')
            #value = float(n[j])  
            value = float(n[1])  #mutation rate
            if value < 1:
                midd[locat].append(value)
    

for i in midd:
    print(i[:10])
#np.random.seed(10)

#data = np.random.randn(8, 4)

print("SO FAR SO GOOD")

new_data = []

for a,b in zip(midd,data):
    new_data.append(a)
    new_data.append(b)

print(new_data[:1])


u_labels = ['1e-7','1e-7 end','1e-6','1e-6 end','1e-5','1e-5 end','1e-4','1e-4 end','0.001','0.001 end','0.005','0.005 end','0.01','0.01 end']

fig, ax = plt.subplots()

ax.set_xticklabels(u_labels,rotation='vertical',fontsize=9)
#ax.set_xticks(u_labels)
ax.set_yscale('log')
ax.set_xlabel('Initial Mutation Rates',fontsize=9)
ax.set_ylabel(opt.ty+' Proliferation Rates',fontsize=9)
i = 0
print(ax.get_xticklabels())
ppl.boxplot(ax, new_data) #,xticklabels=u_labels)

plt.title("Distribution of "+opt.ty+" proliferation rates categorised by initial mutation rate",fontsize=9)
fig.savefig(opt.output+'pro_boxplot.png')
print(opt.output+"DONE")
#ppl.hist(ax,data)
#fig.savefig('histogram_prettyplotlib_default.png')
