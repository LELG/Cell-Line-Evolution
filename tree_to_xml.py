from __future__ import print_function
import random
import os
import math
import sys
import subpopulation

def tree_parse(subpopulation,tumoursize):
    min_size = 1 # 1% min
    F = open(subpopulation.opt["filename"]+"_phylo.xml",'w')
    print("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",file=F)
    print("<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instancei\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">",file=F)
    print("<phylogeny rooted=\"true\">",file=F)
    print("<name>An example</name>",file=F)
    print("<clade>",file=F)
    F.close()
    print_clade(subpopulation, min_size)
    print_end(subpopulation)

def print_end(subpopulation):
    F = open(subpopulation.opt["filename"]+"_phylo.xml",'a')
    print("</clade>",file=F)
    print("</phylogeny>",file=F)
    print("</phyloxml>",file=F)
    F.close()

def print_clade(clade, min_size):  #'clade' is subpop
    if 1: #len(clade.idnt) > 0:
        if clade.size > min_size or clade.depth == 0:
            F = open(clade.opt["filename"]+"_phylo.xml",'a')
            print("<clade branch_length=\"",clade.depth,"\">",file=F)
            print("<name> p:"
                    +str(clade.proliferation)[0:6]+" m:"
                    +str(clade.mutation)[0:6]+" s:"
                    +str(clade.size)+" r:"
                    +str(clade.col)+" "
                    #+clade.idnt
                    +"</name>",file=F)
            F.close()

    for i in range(0,len(clade.nodes)):
        print_clade(clade.nodes[i],min_size)

    if 1: #len(clade.idnt) > 0:
        if clade.size > min_size or clade.depth == 0:
            F = open(clade.opt["filename"]+"_phylo.xml",'a')
            print("</clade>",file=F)
            F.close()


