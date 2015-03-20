from __future__ import print_function
import random
import os
import math
import sys
import subpopulation

def tree_parse(subpop, tumoursize, t_curr, results_dir, fname):
    min_size = 1 # 1% min

    fpath = "{0}/{1}_phylo.xml".format(results_dir, fname)

    xml_header = '''<?xml version="1.0" encoding="UTF-8"?>\n'''
    xml_header += '''<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instancei" '''
    xml_header += '''xsi:schemaLocation="http://www.phyloxml.org '''
    xml_header += '''http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">\n'''
    xml_header += '''<phylogeny rooted="true">\n'''
    xml_header += '''<name>Phylogenetic Tree</name>\n'''
    xml_header += '''<clade>\n'''
    phylo_file = open(fpath, 'w')
    phylo_file.write(xml_header)
    phylo_file.close()
    print_clade(subpop, min_size, float(t_curr), fpath)
    print_xml_footer(fpath)

def print_xml_footer(xml_filepath):
    xml_footer = "</clade>\n"
    xml_footer += "</phylogeny>\n"
    xml_footer += "</phyloxml>"
    phylo_file = open(xml_filepath, 'a')
    phylo_file.write(xml_footer)
    phylo_file.close()

def print_clade(clade, min_size, t_curr, phylo_filepath):  #'clade' is subpop
    if 1: #len(clade.idnt) > 0:
        if clade.size > min_size or clade.depth == 0:
            phylo_file = open(phylo_filepath, 'a')
            if t_curr == 0:
                t_curr = 1
            clade_hdr = '<clade branch_length="{0}">\n'.format(clade.branch_length/t_curr)
            clade_name = "<name> p:{} m:{} s:{} r:{} </name>\n"
            clade_name = clade_name.format(str(clade.prolif_rate)[0:6],
                                           str(clade.mut_rate)[0:6],
                                           str(clade.size),
                                           str(clade.col))
            phylo_file.write(clade_hdr)
            phylo_file.write(clade_name)
            phylo_file.close()

    for subclade in clade.nodes:
        print_clade(subclade, min_size, t_curr, phylo_filepath)

    if 1: #len(clade.idnt) > 0:
        if clade.size > min_size or clade.depth == 0:
            phylo_file = open(phylo_filepath, 'a')
            phylo_file.write("</clade>\n")
            phylo_file.close()
