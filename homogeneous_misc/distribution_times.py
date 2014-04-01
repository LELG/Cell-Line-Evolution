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

for j in range(3):

    lines = [line.strip() for line in open(filename)]
    locat = 0
    for l in lines:
        if l == '' or l == ' ':
            locat+=1
            print("next set")
        else:
            n = l.split(' ')
            print(n)
            #value = float(n[j])  #time
            value = int(n[j])  #time
            data[locat].append(value)
        

    """
    print("data")
    for i in data:
        print(i[:10])
    """

    #np.random.seed(10)

    #data = np.random.randn(8, 4)

    u_labels = ['1e-7','1e-6','1e-5','1e-4','0.001','0.005','0.01','0.05','0.1']

    fig, ax = plt.subplots()

    ax.set_xticklabels(u_labels)
    ax.set_xticks(u_labels)
    ax.set_xlabel('Initial Mutation Rates',fontsize=9)
    ax.set_ylabel('End Time',fontsize=9)
    ax.set_yscale('log')
    i = 0
    print(ax.get_xticklabels())
    ppl.boxplot(ax, data) #,xticklabels=u_labels)

    t = ["time taken to reach full size","time taken until crash","time taken until recovery"]
    tt = ["pre","crash","post"]

    plt.title("distribution of "+t[j]+" categorised by mutation rate",fontsize=9)
    fig.savefig(opt.output+tt[j]+'_boxplot.png')
    print(opt.output+"DONE")
    #ppl.hist(ax,data)
    #fig.savefig('histogram_prettyplotlib_default.png')
