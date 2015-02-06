from __future__ import print_function
import random
import os
import math
import sys
import subpopulation

def tree_parse(subpopulation,tumoursize,time,fname):
    min_size = 1 # 1% min
    F = open(subpopulation.opt["run_dir"]+fname+"_phylo.xml",'w')
    print("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",file=F)
    print("<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instancei\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">",file=F)
    print("<phylogeny rooted=\"true\">",file=F)
    print("<name>Phylogenetic Tree</name>",file=F)
    print("<clade>",file=F)
    F.close()
    print_clade(subpopulation, min_size,float(time),fname)
    print_end(subpopulation,fname)

def print_end(subpopulation,fname):
    F = open(subpopulation.opt["run_dir"]+fname+"_phylo.xml",'a')
    print("</clade>",file=F)
    print("</phylogeny>",file=F)
    print("</phyloxml>",file=F)
    F.close()

def print_clade(clade, min_size,time,fname):  #'clade' is subpop
    if 1: #len(clade.idnt) > 0:
        if clade.size > min_size or clade.depth == 0:
            F = open(clade.opt["run_dir"]+fname+"_phylo.xml",'a')
            print("<clade branch_length=\"",clade.branch_length/time,"\">",file=F)
            print("<name> p:"
                    +str(clade.proliferation)[0:6]+" m:"
                    +str(clade.mutation)[0:6]+" s:"
                    +str(clade.size)+" r:"
                    +str(clade.col)+" "
                    #+clade.idnt
                    +"</name>",file=F)
            F.close()

    for i in range(0,len(clade.nodes)):
        print_clade(clade.nodes[i],min_size,time,fname)

    if 1: #len(clade.idnt) > 0:
        if clade.size > min_size or clade.depth == 0:
            F = open(clade.opt["run_dir"]+fname+"_phylo.xml",'a')
            print("</clade>",file=F)
            F.close()


