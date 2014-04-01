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

labels = ['LightGray','Gray','LightSkyBlue','RoyalBlue','LightSeaGreen','MediumSeaGreen','Khaki','Goldenrod','Tomato','PaleVioletRed','DarkViolet']
data = [[],[],[],[],[],[],[],[],[],[],[]]

lines = [line.strip() for line in open(filename)]
for l in lines:
    n = l.split(' ')
    value = float(n[0]) 
    count = int(n[1])
    group = n[2]
    locat = labels.index(group)
    for x in range(count):
#        if value > 1:
#            value = 1
        if value < 1:
            data[locat].append(value)

#print(data)

#np.random.seed(10)

#data = np.random.randn(8, 4)



u_labels = ['1e-7','1e-6','1e-5','1e-4','0.001','0.005','0.01','0.05','0.1','0.2','0.5']

fig, ax = plt.subplots()

ax.set_xticklabels(u_labels)
ax.set_xticks(u_labels)
ax.set_xlabel('Initial Mutation Rates',fontsize=9)
ax.set_yscale('log')
ax.set_ylabel('End Mutation Rates',fontsize=9)
i = 0
print(ax.get_xticklabels())
for tl in ax.get_xticklabels():
    #tl.set_fontsize(9)
    #tl.set_color(labels[i])
    print(labels[i])
    i+=1
ppl.boxplot(ax, data, xticklabels=u_labels)
fig.savefig(opt.output+'_boxplot.png')
print(opt.output+"DONE")
#ppl.hist(ax,data)
#fig.savefig('histogram_prettyplotlib_default.png')
