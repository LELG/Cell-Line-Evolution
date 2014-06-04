from __future__ import print_function
import os
import time
import subpopulation
import population
import argparse


def drop(subpop, tumoursize, time, label):
    ## print info for multiple populations
    ## print by colour of population
    file_drop = open(subpop.opt["testname"]+"/"+label+"dropdata.dat",'a')
    #summary for all runs - add value and overwrite

    x = subpop.tree_to_list("size_by_col")

    x.sort()

    colours = [0,0,0,0,0,0,0,0,0,0,0]

    #merge colours, super lazynaive method
    #will eventually change whole color sys to be dynamic

    for item in x:
        if item[0] == 'LightGray':
            colours[0] += item[1]
        if item[0] == 'Gray':
            colours[1] += item[1]
        if item[0] == 'LightSkyBlue':
            colours[2] += item[1]
        if item[0] == 'RoyalBlue':
            colours[3] += item[1]
        if item[0] == 'LightSeaGreen':
            colours[4] += item[1]
        if item[0] == 'MediumSeaGreen':
            colours[5] += item[1]
        if item[0] == 'Khaki':
            colours[6] += item[1]
        if item[0] == 'Goldenrod':
            colours[7] += item[1]
        if item[0] == 'Tomato':
            colours[8] += item[1]
        if item[0] == 'PaleVioletRed':
            colours[9] += item[1]
        if item[0] == 'DarkViolet':
            colours[10] += item[1]


    output = ""
    for i in colours:
        output = output + " " + str(i)
    print(output, file = file_drop)

def read_drop(fn):
    file_mid = open(fn+"/sum_hetpop_middist.dat",'a')
    file_end = open(fn+"/sum_hetpop_enddist.dat",'a')
    #MID
    print("reading...")
    colours = [0,0,0,0,0,0,0,0,0,0,0]
    #Read the info for all colour populations, and graph
    filename = fn+"/"+"mid"+"dropdata.dat"
    lines = [line.strip() for line in open(filename)]
    print(lines)
    #sum total counts
    for l in lines:
        aa = l.split(' ')
        c = [int(z) for z in aa] 
        print("and then")
        colours = [x + y for x, y in zip(colours, c)]
    print("ALL MID COLOURS ",colours, file = file_mid)
    #then plot it

    #END
    print("reading...")
    colours = [0,0,0,0,0,0,0,0,0,0,0]
    filename = fn+"/"+"end"+"dropdata.dat"
    lines = [line.strip() for line in open(filename)]
    #sum total counts
    for l in lines:
        aa = l.split(' ')
        c = [int(z) for z in aa] 
        colours = [x + y for x, y in zip(colours, c)]
    print("ALL END COLOURS ",colours, file = file_end)
    #then plot it
    

     
#parser = argparse.ArgumentParser()
#parser.add_argument('-f','--filez', default='filez')
#opt = parser.parse_args()
#read_drop(opt.filename,'mid')
#read_drop(opt.filename,'end')


