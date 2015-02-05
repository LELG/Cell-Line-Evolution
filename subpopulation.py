from __future__ import print_function
import random
import os
import csv
import math
import sys
import tree_to_xml
import numpy.random
#import tree_to_xml_screen

class Subpopulation():
    def __init__(self,opt,p,m,depth,time,mut_type,col,pmp,pmn,pim,pdm,msc,prev_time):
        self.opt = opt #convert opt to dictionary
        self.proliferation = p
        self.mutation = m
        self.mut_static = m #static for scale calution of new value
        self.mscale = msc
        self.death = opt["die"]
        self.size = 1
        self.precrash_size = 0
        self.maxsize_lim = opt["maxsize_lim"]
        self.prolif_lim = opt["prolif_lim"]
        self.prolif_adj = opt["prolif_lim"]
        self.nodes = []
        self.depth = depth
        self.s_time = time
        self.d_time = self.opt["max_cycles"]+1
        self.idnt = "" ##This gets created when we run freq through tree
        self.prob_mut_pos = pmp #opt["prob_mut_pos"]
        self.prob_mut_neg = pmn #opt["prob_mut_neg"]
        self.mut_type = mut_type  #b/d/n
        self.prob_inc_mut = pim #opt["prob_inc_mut"]
        self.prob_dec_mut = pdm #opt["prob_dec_mut"]
        self.mutagenic_pressure = 0
        self.col = col
        self.branch_length = time - prev_time

    def newsubpop_from_file(self, filename):
        """Load a heterogeneous starting population from a subfile."""
        subfile = open(filename)
        reader = csv.DictReader(subfile)
        for line in reader:
            mut = float(reader["#mut"])
            pmp = float(reader["pm+"])
            pmn = float(reader["pm-"])
            pim = float(reader["pim"])
            pdm = float(reader["pdm"])
            msc = float(reader["msc"])
            col = reader["col"]

            #add new subpopulation node
            new_subpop = Subpopulation(self.opt, \
                    self.proliferation, mut,
                    1, 0, 'n', col, pmp, pmn, pim, pdm, msc, 0)
            new_subpop.size = self.opt["init_size"]
            self.nodes.append(new_subpop)

    def info(self):
        print("Population Info")
        print(self.opt)
        print("No. of nodes: ",len(self.nodes))
        print(self.proliferation,self.mutation)
        for i,j in enumerate(self.nodes):
            self.nodes[i].info()

    def set_precrash_size(self):
        """ Record pre crash size """
        self.precrash_size = self.size
        for elem in self.nodes:
            elem.set_precrash_size()

    def cycle(self, tumoursize, select_pressure, mutagenic_pressure, time, prolif_adj):
        self.mutagenic_pressure = mutagenic_pressure
        cells_mut = 0
        if self.size <= 0: #DEAD
            self.size = 0
            if self.d_time == self.opt["max_cycles"]+1: #replace with real time of death
                self.d_time = time

        else:
            if self.opt["Z"]: #if running large dataset
            # MIGHT BE QUICKER TO FILTER ALL END NODES
            # ONE CALL PER 10 RUNS / 2 RUNS 
            # remove 0 nodes
            #MERGE ALIVE INTO CYCLE
            #EATS UP THE ACTUAL REFERENCE FOR TREE
            #CHOMPS ANY INTERMEDIATE NODES
                new_nodes_1 = []
                for i in range(0, len(self.nodes)): #iterate over copy of list
                    #keep if > 0 
                    if self.nodes[i].size > 0:
                        new_nodes_1 = new_nodes_1 + [self.nodes[i]]
                    else:
                        #keep if it has subsequent nodes
                        if self.nodes[i].nodes:
                            new_nodes_1 = new_nodes_1 + [self.nodes[i]]

                self.nodes = new_nodes_1
        
            #IF ALIVE
            self.prolif_adj = prolif_adj
            new_p = self.proliferation - self.prolif_adj - select_pressure
            if self.mutagenic_pressure != 0:
                new_m = self.mutation * self.mutagenic_pressure
            else:
                new_m = self.mutation
            cells_dead = numpy.random.binomial(self.size, self.death)
            cells_new = 0
            if new_p > 0:
                cells_new = numpy.random.binomial(self.size, new_p) #IF THIS SMALL?? shouldnt get neg
            #else:
            #    print("IT HAPPENDED AGAIN")
            cells_mut = 0
            if cells_new > 0 and new_m > 0:
                cells_mut = numpy.random.binomial(cells_new, new_m)

            self.size = self.size + cells_new - cells_mut - cells_dead

        """ Repeat for subsequent nodes """
        new_sizes = cells_mut
        sub_count = 0
        if self.size > 0:
            #mut_avg = ((self.mutation + mutagenic_pressure) * self.size)
            #check this
            mut_avg = new_m * self.size
            pro_avg = ((self.proliferation-self.prolif_adj) * self.size)
        else:
            mut_avg = 0
            pro_avg = 0
        for i in range(0, len(self.nodes)):
            ret_new_size, ret_sub_count, ret_mut_avg, ret_pro_avg \
                                        = self.nodes[i].cycle(tumoursize, \
                                          select_pressure, mutagenic_pressure,
                                          time, prolif_adj)
            #COUNT NEW CELLS
            new_sizes = new_sizes + ret_new_size
            #COUNT NO. CLONES
            sub_count = sub_count + ret_sub_count
            #COUNT MUT AVG
            mut_avg = mut_avg + ret_mut_avg
            #COUNT PRO AVG
            pro_avg = pro_avg + ret_pro_avg

            #COUNT PRO AVG

        if cells_mut > 0: #make new subpopulations AFTER THE FACT 
            for i in range(0, cells_mut):
                self.new_subpop(time)

        if self.size > 0:
            sub_count += 1
        #also return 
        return self.size + new_sizes, sub_count, mut_avg, pro_avg
        #return new tumour size, new number of clones


    """
    def mutate(self, time, mutagenic_pressure):
        y = random.random()
        if y < float(self.mutation):
            self.new_subpop(time)
            return 1
        return 0
    """
            
    def new_subpop(self, time):
        #AMOUNT OF CHANGE
        #print("mm")
        mut_type = 'n'
        x = random.random()
        p = self.proliferation
        if x < self.prob_mut_pos: #NEG
            p = self.mutation_beneficial()
            mut_type = 'b'

        if (self.prob_mut_pos + self.prob_mut_neg) > x > self.prob_mut_pos:
            p = self.mutation_deleterious()
            mut_type = 'd'

        z = random.random()
        m = self.mutation
        if z < self.prob_inc_mut:
            m = self.mutation_change_increase()
        if (self.prob_inc_mut + self.prob_dec_mut) > z > self.prob_inc_mut:
            m = self.mutation_change_decrease()

        size = 1
        depth = self.depth + 1
        #new_opt = self.opt
        col = self.col
        pmp = self.prob_mut_pos
        pmn = self.prob_mut_neg
        pim = self.prob_inc_mut
        pdm = self.prob_dec_mut
        msc = self.mscale
        self.nodes.append(Subpopulation(self.opt, p, m, depth, time, mut_type, col, \
                pmp, pmn, pim, pdm, msc, self.s_time)) #add new subpop

    """ GATHER DATA FOR ANALYTICS """
    """ Dangerous loops beyond this sign """

    def mutation_beneficial(self):
        """ Gives random value between 0-1

        Scale the value down to a max of 0.001
        OR dynamically the max could be difference between P and D
        OR one unit of 'selective pressure'

        """
        alpha = 1
        beta = 3
        SCALE = self.opt["scale"] * self.opt["pro"] 
        new_proliferation_delta =\
                numpy.random.beta(alpha,beta,size=1)[0]*SCALE
        p = self.proliferation + new_proliferation_delta
        #print("CHANGE PRO ",self.proliferation," to ",p, " using ", new_proliferation_delta, " with scale ",SCALE)
        if p > 1:
            p = 0.99
        return p

    def mutation_deleterious(self):
        """ Gives random value between 0-1

        Scale the value down to a max of 0.001
        OR dynamically the max could be difference between P and D
        OR one unit of 'selective pressure'

        """
        alpha = 1
        beta = 3
        SCALE = self.opt["scale"] * self.opt["pro"] 
        new_proliferation_delta =\
                numpy.random.beta(alpha,beta,size=1)[0]*SCALE
        p = self.proliferation - new_proliferation_delta
        if p < 0:
            p = SCALE / 10000.0
        return p
        
    def mutation_change_increase(self):
        alpha = 1
        beta = 3
        #SCALE = 10000 
        #SCALE = 100000000
        #print("SCALE",self.opt["mscale"]," and ",self.opt["mut"])
        SCALE = self.mscale * self.mutation
        new_mutation_delta = numpy.random.beta(alpha,beta,size=1)[0]*SCALE
        m = self.mutation + new_mutation_delta
        #print("CHANGE ",self.mutation," to ",m, " using ", new_mutation_delta, " with scale ",SCALE)
        if m > 1:
            m = 0.99
        return m
        
    def mutation_change_decrease(self):
        alpha = 1
        beta = 3
        #SCALE = 10000 
        #SCALE = 100000000
        SCALE = self.mscale * self.mutation
        new_mutation_delta = numpy.random.beta(alpha,beta,size=1)[0]*SCALE
        m = self.mutation - new_mutation_delta
        if m < 0:
            m = SCALE / 10000.0
        return m
        
    def pop_as_list(self):
        n = self.tree_to_list("size")
        n.sort()
        #out = filter(lambda a: a!=0, n)
        return n

    def pop_count(self):
        x = self.pop_as_list()
        z = 0 
        for i in x:
            z = z + i
        return z

    def freq_of_mutation(self,tumoursize):
        x = self.freq_to_list("0-0")
        y = []
        z = []
        for i, j in x:
            y = y + [i/float(tumoursize), j]
            z = z + [i/float(tumoursize)]
        return x, y, z

    def freq_to_list(self,idnt):
        self_node = []
        blank_ctree = []
        if self.size > 0:
            self_node = [(self.pop_count(), \
                    "pr-"+str(self.proliferation)+"-"+ \
                    self.mut_type+idnt)]
            self.idnt = idnt + str(self.depth)
            #if len(self.nodes) == 0: #deepest point in tree
                #self_node = [
        
        if self.nodes == []:
            return self_node

        else:
            for i in range(0,len(self.nodes)):
                blank_ctree = blank_ctree + \
                        self.nodes[i].freq_to_list(str(i)+self.nodes[i].mut_type+idnt)
        return blank_ctree + self_node

    """ Instead of one function for each kind, better to return """
    """ list of objects, right? """

    def tree_to_list(self,prm): 
        """ SIZE AS FILTER """
        #make it get everything, then use dictionary to access field?
        #TO DO have size optional paramter to filter with 
        self_node = []
        blank_ctree = []
        if prm == "size":
            if self.size > 0:
                self_node = [self.size]

        if prm == "size_by_col":
            if self.size > 0:
                self_node = [(self.col, self.size)]

        if prm == "proliferation":
            if self.size > 0:
                self_node = [self.proliferation]

        if prm == "proliferation_size":
            if self.size > 0:
                self_node = [(self.proliferation, self.size)]

        if prm == "mutation":
            if self.size > 0:
                if self.mutagenic_pressure > 0:
                    self_node = [self.mutation*self.mutagenic_pressure] 
                else:
                    self_node = [self.mutation]
                    #added mutagenic pressure

        if prm == "effective_proliferation":
            if self.size > 0:
                self_node = [(self.proliferation-self.prolif_adj,self.size)] 

        #could change to effective mutation rate for mutagenic selection
        if prm == "mutation_rate":
            if self.size > 0:
                self_node = [(self.mutation,self.size)]

        if prm == "cell_line_time":
            self_node = [(self.col,self.s_time,self.d_time)]
                #select pr?

        if prm == "cell_line_time_mut":
            self_node = [(self.s_time,self.d_time,self.mutation,self.size)]
                #select pr?

        if prm == "circles":
            if self.size > 0:
                self_node = [(self.mutation,self.proliferation,self.size)]

        if prm == "circles_all":
            self_node = [(self.mutation,self.proliferation,self.size)]

        if prm == "mutation_distribution":
            if self.size > 0:
                self_node = [(self.mutation,self.precrash_size,self.size)]

        if prm == "two_side_size":
            if self.size > 0:
                self_node = [(self.precrash_size,self.size)]

        if prm == "mutation_distribution_1":
            if self.size > 0:
                self_node = [(self.mutation,self.proliferation,self.size,self.precrash_size)]

        if self.nodes == []:
            return self_node
        else:
            for i in range(0,len(self.nodes)):
                 blank_ctree = blank_ctree + self.nodes[i].tree_to_list(prm)
        return blank_ctree + self_node
