"""
Module for classes and functions related to
individual clones.

Authors
-------
Andrew Bakshi : andrew.bakshi@gmail.com
Yoshua Wakeham : yoshwakeham@gmail.com
"""
from __future__ import print_function
import random
import csv
import json
import numpy as np

class Subpopulation(object):
    """
    Individual cancer clone.
    """
    def __init__(self, opt, p, m, depth, time, mut_type, col, pmp, pmn, pim, pdm, prev_time):
        """Create new clone."""
        self.opt = opt
        self.proliferation = p
        self.mutation = m
        self.mut_static = m #static for scale calution of new value
        self.death = opt.die
        self.size = 1
        self.precrash_size = 0
        self.max_size_lim = opt.max_size_lim
        # initialise prolif_adj to maximum amount
        self.prolif_adj = opt.prolif_lim
        self.nodes = []
        self.depth = depth
        self.s_time = time
        self.d_time = None
        self.idnt = "" ##This gets created when we run freq through tree
        self.prob_mut_pos = pmp #opt["prob_mut_pos"]
        self.prob_mut_neg = pmn #opt["prob_mut_neg"]
        self.mut_type = mut_type  #b/d/n
        self.prob_inc_mut = pim #opt["prob_inc_mut"]
        self.prob_dec_mut = pdm #opt["prob_dec_mut"]
        self.mutagenic_pressure = 0
        self.col = col
        self.branch_length = time - prev_time

    def newsubpop_from_file(self, filename, init_size):
        """Load a heterogeneous starting population from a subfile."""
        subfile = open(filename)
        reader = csv.DictReader(subfile)
        for line in reader:
            mut = float(line["#mut"])
            pmp = float(line["pm+"])
            pmn = float(line["pm-"])
            pim = float(line["pim"])
            pdm = float(line["pdm"])
            msc = float(line["msc"])
            col = line["col"]

            #add new subpopulation node
            new_subpop = Subpopulation(opt=self.opt,
                                       p=self.proliferation,
                                       m=mut, depth=1, time=0,
                                       mut_type='n', col=col,
                                       pmp=pmp, pmn=pmn,
                                       pim=pim, pdm=pdm,
                                       prev_time=0)
            new_subpop.size = init_size
            self.nodes.append(new_subpop)

    def info(self):
        """Print info about this clone and all its children"""
        print("Clone Info")
        print(self.opt)
        print("No. of nodes: ", len(self.nodes))
        print(self.proliferation, self.mutation)
        for node in self.nodes:
            node.info()

    def set_precrash_size(self):
        """ Record pre crash size """
        self.precrash_size = self.size
        for node in self.nodes:
            node.set_precrash_size()

    def update(self, select_pressure, mutagenic_pressure, time, prolif_adj, prune_clones=False):
        """Update this clone and its children for one time step."""
        self.mutagenic_pressure = mutagenic_pressure
        self.prolif_adj = prolif_adj
        cells_new = 0
        cells_mut = 0

        if self.is_dead():
            self.size = 0
            if not self.d_time:
                self.d_time = time
        else:
            if prune_clones: #if running large dataset
            # MIGHT BE QUICKER TO FILTER ALL END NODES
            # ONE CALL PER 10 RUNS / 2 RUNS
            # remove 0 nodes
            #MERGE ALIVE INTO CYCLE
            #EATS UP THE ACTUAL REFERENCE FOR TREE
            #CHOMPS ANY INTERMEDIATE NODES
                nodes_to_keep = [node for node
                                 in self.nodes
                                 if not node.is_dead_end()]
                self.nodes = nodes_to_keep

            effective_pro = self.proliferation - prolif_adj - select_pressure

            effective_mut = self.mutation
            if self.mutagenic_pressure:
                effective_mut *= self.mutagenic_pressure

            # sample for number of dead cells
            cells_dead = np.random.binomial(self.size, self.death)
            # sample for new cells
            if effective_pro > 0:
                cells_new = np.random.binomial(self.size, effective_pro)
            # sample for mutations
            if cells_new > 0 and effective_mut > 0:
                cells_mut = np.random.binomial(cells_new, effective_mut)

            self.size += cells_new - cells_mut - cells_dead

        # update children

        new_pop_size = cells_mut
        sub_count = 0
        if not self.is_dead():
            mut_agg = effective_mut * self.size
            pro_agg = ((self.proliferation - self.prolif_adj) * self.size)
        else:
            mut_agg = 0
            pro_agg = 0

        for node in self.nodes:
            node_results = node.update(select_pressure, mutagenic_pressure,
                                       time, prolif_adj)
            ret_new_size, ret_sub_count, ret_mut_agg, ret_pro_agg = node_results

            new_pop_size += ret_new_size
            sub_count += ret_sub_count
            mut_agg += ret_mut_agg
            pro_agg += ret_pro_agg

        # make new subpopulations after the fact
        for _ in range(0, cells_mut):
            self.new_subpop(time)

        if self.size > 0:
            sub_count += 1
            new_pop_size += self.size

        return new_pop_size, sub_count, mut_agg, pro_agg

    def new_subpop(self, time):
        """Spawn a new child clone."""
        # get child's proliferation rate
        prolif_mut_event = random.random()

        if prolif_mut_event <= self.prob_mut_pos:
            new_prolif = self.mutation_beneficial()
            new_mut_type = 'b'
        elif prolif_mut_event <= self.prob_mut_pos + self.prob_mut_neg:
            new_prolif = self.mutation_deleterious()
            new_mut_type = 'd'
        else:
            new_prolif = self.proliferation
            new_mut_type = 'n'

        # get child's mutation rate
        mut_mut_event = random.random()

        if mut_mut_event <= self.prob_inc_mut:
            new_mut = self.mutation_change_increase()
        elif mut_mut_event <= self.prob_inc_mut + self.prob_dec_mut:
            new_mut = self.mutation_change_decrease()
        else:
            new_mut = self.mutation

        size = 1
        new_depth = self.depth + 1
        new_child = Subpopulation(opt=self.opt,
                                  p=new_prolif, m=new_mut,
                                  depth=new_depth, time=time,
                                  mut_type=new_mut_type, col=self.col,
                                  pmp=self.prob_mut_pos, pmn=self.prob_mut_neg,
                                  pim=self.prob_inc_mut, pdm=self.prob_dec_mut,
                                  prev_time=self.s_time)
        self.nodes.append(new_child)

    def mutation_beneficial(self):
        """Gives random value between 0-1

        Scale the value down to a max of 0.001
        OR dynamically the max could be difference between P and D
        OR one unit of 'selective pressure'
        """
        alpha = 1
        beta = 3
        pro_scale = self.opt.scale * self.opt.pro
        prolif_delta = np.random.beta(alpha, beta, size=1)[0] * pro_scale
        new_prolif = self.proliferation + prolif_delta
        if new_prolif > 1:
            new_prolif = 0.99
        return new_prolif

    def mutation_deleterious(self):
        """ Gives random value between 0-1

        Scale the value down to a max of 0.001
        OR dynamically the max could be difference between P and D
        OR one unit of 'selective pressure'
        """
        alpha = 1
        beta = 3
        pro_scale = self.opt.scale * self.opt.pro
        prolif_delta = np.random.beta(alpha, beta, size=1)[0] * pro_scale
        new_prolif = self.proliferation - prolif_delta
        if new_prolif < 0:
            new_prolif = pro_scale / 10000.0
        return new_prolif

    def mutation_change_increase(self):
        alpha = 1
        beta = 3
        mut_scale = self.opt.mscale * self.mutation
        mut_delta = np.random.beta(alpha, beta, size=1)[0] * mut_scale
        new_mut = self.mutation + mut_delta
        if new_mut > 1:
            new_mut = 0.99
        return new_mut

    def mutation_change_decrease(self):
        alpha = 1
        beta = 3
        mut_scale = self.opt.mscale * self.mutation
        mut_delta = np.random.beta(alpha, beta, size=1)[0] * mut_scale
        new_mut = self.mutation - mut_delta
        if new_mut < 0:
            new_mut = mut_scale / 10000.0
        return new_mut

    # GATHER DATA FOR ANALYTICS
    # Dangerous loops beyond this sign

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

    def freq_to_list(self, idnt):
        self_node = []
        blank_ctree = []
        if self.size > 0:
            self_node = [(self.pop_count(),
                          "pr-{}-{}{}".format(str(self.proliferation),
                                              self.mut_type, idnt))]
            self.idnt = idnt + str(self.depth)
            #if len(self.nodes) == 0: #deepest point in tree
                #self_node = [

        if not self.nodes:
            return self_node
        else:
            for node in self.nodes:
                new_idnt = str(self.nodes.index(node)) + node.mut_type + idnt
                blank_ctree += node.freq_to_list(new_idnt)
            return blank_ctree + self_node

    # Instead of one function for each kind, better to return
    # list of objects, right?

    def tree_to_list(self, prm):
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
                    self_node = [self.mutation * self.mutagenic_pressure]
                else:
                    self_node = [self.mutation]
                    #added mutagenic pressure

        if prm == "effective_proliferation":
            if self.size > 0:
                self_node = [(self.proliferation - self.prolif_adj, self.size)]

        #could change to effective mutation rate for mutagenic selection
        if prm == "mutation_rate":
            if self.size > 0:
                self_node = [(self.mutation, self.size)]

        if prm == "cell_line_time":
            self_node = [(self.col, self.s_time, self.d_time)]

        if prm == "cell_line_time_mut":
            self_node = [(self.s_time, self.d_time, self.mutation, self.size)]

        if prm == "circles":
            if self.size > 0:
                self_node = [(self.mutation, self.proliferation, self.size)]

        if prm == "circles_all":
            self_node = [(self.mutation, self.proliferation, self.size)]

        if prm == "mutation_distribution":
            if self.size > 0:
                self_node = [(self.mutation, self.precrash_size, self.size)]

        if prm == "two_side_size":
            if self.size > 0:
                self_node = [(self.precrash_size, self.size)]

        if prm == "mutation_distribution_1":
            if self.size > 0:
                self_node = [(self.mutation, self.proliferation,
                              self.size, self.precrash_size)]

        if not self.nodes:
            return self_node
        else:
            for node in self.nodes:
                blank_ctree += node.tree_to_list(prm)
            return blank_ctree + self_node

    def is_dead(self):
        return self.size <= 0

    def has_children(self):
        return bool(self.nodes)

    def is_dead_end(self):
        return self.is_dead() and not self.has_children()

    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,
                          sort_keys=True, indent=4)
