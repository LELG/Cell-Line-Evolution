from __future__ import print_function
import numpy
import argparse
import prettyplotlib as ppl
from prettyplotlib import plt


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', default='filename')
opt = parser.parse_args()

filename = opt.filename

"""
fig = pyplot.figure()
ax = fig.add_subplot(111)

N = 7
ind = numpy.arange(N)
width = 0.5

#change these according to results
#BASE_0

succ_runs = 42.0
vals = [0.001, 0.001, 1, 1, 5, 21, 15]
print(vals)
vals = [(i / succ_runs)*100 for i in vals]
print(vals)
colors = ['k','c','b','m','r','orange','chartreuse']

ax.bar(ind, vals, width, color=colors)
pyplot.xticks(ind+width/2.,('0.001','0.005','0.001','0.005','0.01','0.05','0.1'))
pyplot.ylabel("Recovery Percent")
pyplot.title("Composition of Populations which Developed Resistance")
pyplot.savefig("bar_mut")
pyplot.clf()

"""

###
#FORCE RANGE TO 0 - 100

fig, ax = plt.subplots(1) 

N = 11
ind = numpy.arange(N)
print(ind)
width = 0.5

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}

#          k      c      b      m       r       o        char
vals = [ 45174 , 48810 , 41719 , 255017 , 642090 , 3053979 , 3066135 , 1731655 , 582226 , 126045 , 6738 ]
total = float(sum(vals))
print(vals)
vals = [(i / total)*100 for i in vals]
print(vals)
colors = ['LightGray','Gray','LightSkyBlue','RoyalBlue','LightSeaGreen','MediumSeaGreen','Khaki',\
        'Goldenrod','Tomato','MediumVioletRed','DarkViolet']

#plt.bar(ind, vals, width)
ppl.bar(ax, ind, vals, color=colors)
plt.xticks(ind+width/2.,\
        ('1e-7','1e-6','1e-5','1e-4','0.001','0.005','0.01', '0.05', '0.1', '0.2', '0.5')\
        ,fontsize=11)
plt.yticks(fontsize=11)
plt.title("Distribution of cells by mutation rate across resistant populations",fontsize=11)
#plt.ylim((0,100))
plt.ylabel("Percent of Total",fontsize=11)
plt.xlabel("Initial Mutation rate: u",fontsize=11)
plt.savefig(opt.filename+"_bar")
plt.clf()


###
