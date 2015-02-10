"""
High level simulation routines.

Authors
-------
Andrew Bakshi : andrew.bakshi@gmail.com
Yoshua Wakeham : yoshua.wakeham@petermac.org
"""

from __future__ import print_function
import os
import csv
from analytics import Analytics
from subpopulation import Subpopulation
import tree_to_xml
import time
import dropdata
from utilities import secs_to_hms


class Population(object):
    """
    Contains high level simulation routines.

    Create initial clone (subpopulation) from input parameters
    Run simulation
    Store data in analytics for each loop
    Print summary data to file
    Print graphs
    """
    def __init__(self, opt):
        self.opt = opt
        self.tumoursize = opt.init_size
        self.clonecount = 1
        self.max_size_lim = opt.max_size_lim
        self.prolif_lim = self.opt.prolif_lim
        depth = 0
        time = 0
        mut_type = 'n'
        self.subpop = Subpopulation(vars(self.opt),
                                    self.opt.pro, self.opt.mut,
                                    depth, time, mut_type, 'n',
                                    self.opt.prob_mut_pos,
                                    self.opt.prob_mut_neg,
                                    self.opt.prob_inc_mut,
                                    self.opt.prob_dec_mut,
                                    self.opt.mscale, 0)
        self.subpop.size = opt.init_size
        self.analytics_base = Analytics()
        self.select_pressure = 0.0
        self.mutagenic_pressure = 0.0
        self.selective_pressure_applied = False
        self.avg_pro_rate = self.opt.pro
        self.avg_mut_rate = self.opt.mut
        self.mid_proliferation = []
        self.mid_mutation = []
        if opt.init_diversity:
            self.subpop.size = 0 #don't use init dummy popn if reading from file
            self.subpop.newsubpop_from_file(self.opt.sub_file)

    def info(self):
        print("Parameter Set: ", self.opt)

    def update(self, treatmt, t):
        """
        Recaculate ratio of population size
        Call analytics
        Call subpopulation cycle for each subpopulation
        (Calculate 1 discrete time step for each subpopulation)
        """
        if treatmt.is_introduced:
            # update amount of selective and mutagenic pressure
            self.select_pressure = treatmt.select_pressure
            self.mutagenic_pressure = treatmt.mutagenic_pressure

        # determine proliferation adjustment due to
        # population size
        ratio = self.tumoursize / float(self.max_size_lim)
        prolif_adj = ratio * float(self.prolif_lim)

        # update subpopulations, getting back
        # tumour size, clone count, and aggregated
        # mutation and proliferation rates
        subpop_results = self.subpop.cycle(self.tumoursize,
                                           self.select_pressure,
                                           self.mutagenic_pressure,
                                           t, prolif_adj)
        self.tumoursize, self.clonecount, agg_mut, agg_pro = subpop_results

        if not self.is_dead():
            self.avg_mut_rate = agg_mut / float(self.tumoursize)
            self.avg_pro_rate = agg_pro / float(self.tumoursize)

    def is_dead(self):
        return self.tumoursize <= 0

    def exceeds_size_limit(self, max_size_lim, tolerance=0.0):
        """Determine whether this population has exceeded a certain size."""
        size_limit = max_size_lim * (1 + tolerance)
        return self.tumoursize > size_limit
