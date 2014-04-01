from __future__ import print_function
import numpy
from matplotlib import pyplot
import argparse
from matplotlib import rc
import prettyplotlib as ppl


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

fig = pyplot.figure()
ax = fig.add_subplot(111)

N = 7
ind = numpy.arange(N)
width = 0.5

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}

rc('font', **font)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
pyplot.rc('font', family='serif')

#change these according to results
#
#          k      c      b      m       r       o        char
vals = [ 202916 , 1162586 , 1996777 , 3333848 , 2243760 , 967952 , 854 ]
total = float(sum(vals))
print(vals)
vals = [(i / total)*100 for i in vals]
print(vals)
colors = ['k','c','b','m','r','orange','chartreuse']

ax.bar(ind, vals, width, color=colors)
pyplot.xticks(ind+width/2.,('0.000001','0.000005','0.00001','0.00005','0.0001','0.0005','0.01'))
pyplot.title("Total number of cells across resistant populations")
pyplot.ylim((0,100))
pyplot.ylabel("Composition of Populations which Developed Resistance")
pyplot.xlabel("Mutation rate: u",fontsize=10)
pyplot.savefig(opt.filename+"_bar")
pyplot.clf()


###
