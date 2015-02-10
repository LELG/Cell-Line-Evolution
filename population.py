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

    def cycle(self, t):
        """
        Recaculate ratio of population size
        Call analytics
        Call subpopulation cycle for each subpopulation
        (Calculate 1 discrete time step for each subpopulation)
        """
        ratio = self.tumoursize / float(self.max_size_lim)
        prolif_adj = ratio * float(self.prolif_lim)

        subpop_results = self.subpop.cycle(self.tumoursize,
                                           self.select_pressure,
                                           self.mutagenic_pressure,
                                           t, prolif_adj)
        self.tumoursize, self.clonecount, agg_mut, agg_pro = subpop_results
        if self.tumoursize > 0:
            self.avg_mut_rate = agg_mut / float(self.tumoursize)
            self.avg_pro_rate = agg_pro / float(self.tumoursize)

    def print_results(self, when, end_time):
        """ Print all results to plots / file

        Print result to graphs using matplotlib
        If --R option parsed print raw output to file

        """

        filename = "{0}/{1}/{2}/{3}".format(self.opt.test_group_dir,
                                            self.opt.param_set,
                                            self.opt.run_number,
                                            when)

        if self.opt.R:
            from outputdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat
        else:
            from plotdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat,\
                    mutation_distribution, mutation_crash, \
                    make_subpop_life_mut, make_popsub

        # Print only once at end of simulation #

        if self.selective_pressure_applied: #ONLY PRINT AT END OF SIM
            #POPULATION
            make_plot(self.analytics_base.time,
                      self.analytics_base.population,
                      filename + "population_graph", "Population Size")
            #SUBPOPULATION
            make_plot(self.analytics_base.time,
                      self.analytics_base.subpopulation,
                      filename + "subpop_graph", "No. of Clones")

            #POPULATION + SUBPOPULATION
            make_popsub(self.analytics_base.time,
                        self.analytics_base.population,
                        self.analytics_base.time,
                        self.analytics_base.subpopulation,
                        filename + "popsubpop", 'Tumour Size', 'No. of Clones')
            #MUTATION RATE BUT SUBPOPULATION
            make_popsub(self.analytics_base.time,
                        self.analytics_base.mutation,
                        self.analytics_base.time,
                        self.analytics_base.subpopulation,
                        filename + "mutsubpop", 'Average Mutation Rate',
                        'No. of Clones')

            #PROLIFERATION RATE + MUTATION RATE
            make_popsub(self.analytics_base.time,
                        self.analytics_base.mutation,
                        self.analytics_base.time,
                        self.analytics_base.proliferation,
                        filename + "prolifmut", 'Mutation Rate',
                        'Proliferation Rate')
            #PROLIFERATION RATE + POPULATION
            make_popsub(self.analytics_base.time,
                        self.analytics_base.population,
                        self.analytics_base.time,
                        self.analytics_base.proliferation,
                        filename + "prolifandpop", 'Tumour Size',
                        'Proliferation Rate')

            #EFFECTIVE PROLIFERATION RATE
            make_plot(self.analytics_base.time,
                      self.analytics_base.proliferation,
                      filename + "effect_prolif", "EFFECT PROLIFERATION")
            #MUTATION RATE AVG
            make_plot(self.analytics_base.time,
                      self.analytics_base.mutation,
                      filename + "mutation_avg", "MUTATION AVG")

            # Mutation...
            mut_distro = self.subpop.tree_to_list("mutation_distribution")
            mutation_distribution(mut_distro,
                                  filename + "mutation_distribution",
                                  "MUTATION V TIME - PRE/POST CRASH",
                                  self.opt.scale)

            # Mutation
            if self.tumoursize > 0:
                mut_distro = self.subpop.tree_to_list("mutation_distribution_1")
                mutation_crash(mut_distro,
                               filename + "mutation_distribution_1",
                               "MUTATION V TIME - PRE/POST CRASH",
                               self.opt.scale)

            #cell lines graph  - [(self.s_time,self.d_time)]
            make_subpop_life(self.subpop.tree_to_list("cell_line_time"),
                             filename + "cell_lines_alpha",
                             "CELL LIFESPAN", end_time, self.opt.max_cycles,
                             self.opt.select_time)

            """
            #cell lines graph  - [(self.s_time,self.d_time)]
            make_subpop_life_mut(self.subpop.tree_to_list("cell_line_time_mut"),
                    filename+"cell_lines_muta",
                    "CELL LIFESPAN",end_time,self.opt.max_cycles,
                    self.opt.select_time,self.opt.mut,self.tumoursize)
            """


        # Print at mid and end of simulation #

        #PRINT END HISTOGRAM IF POPULATION STILL ALIVE
        if self.tumoursize > 0:
            #PROLIFERATION HISTOGRAM
            # [(self.proliferation-self.prolif_adj,self.size)]

            pro_hist = self.subpop.tree_to_list("proliferation")
            pro_hist.sort()
            mut_hist = self.subpop.tree_to_list("mutation")
            mut_hist.sort()

            make_hist(pro_hist,
                      filename+"proliferation_hist",
                      "PROLIFERATION RATES", 25)
            #MUTATION HISTOGRAM
            # [self.mutation]
            make_hist(mut_hist,
                      filename+"mutation_hist",
                      "MUTATION RATES", 25)
            #POPULATION HISTOGRAM
            # [self.size]

            pop_hist = self.subpop.tree_to_list("size"),

            make_hist(pop_hist,
                      filename+"population_hist",
                      "POPULATION DIVISION", 25)
            #ALLELE FREQ
            # just_allele_freq z = z + [i/float(tumoursize)]
            norm, r, just_allele_freq = \
                    self.subpop.freq_of_mutation(self.tumoursize)
            make_hist(just_allele_freq,
                      filename + "allele",
                      "ALLELE FREQ",
                      #bins equal to number of sub pops
                      self.analytics_base.subpopulation[-1])
            #CELL CIRCLE MUT V PRO RATES
            # if size > 0 [(self.mutation,self.proliferation,self.size)]
            mutation_v_proliferation(self.subpop.tree_to_list("circles"),
                                     filename + "circles",
                                     "MUTATION V PROLIFERATION RATES",
                                     self.opt.scale)

            # [(self.mutation,self.proliferation,self.size)]
            mutation_v_proliferation(self.subpop.tree_to_list("circles_all"),
                                     filename + "circles_all",
                                     "MUTATION V PROLIFERATION RATES",
                                     self.opt.scale)

            #MAKE CIRCLES ACROSS ALL GRAPHS BY WRITING TO 1 FILE
            # if size > 0 [(self.mutation,self.proliferation,self.size)]
            fpath = "{0}/{1}/{1}-circles.dat".format(self.opt.test_group_dir,
                                                     self.opt.param_set)
            mutation_v_proliferation_dat(self.subpop.tree_to_list("circles"),
                                         fpath,
                                         "MUTATION V PROLIFERATION RATES",
                                         self.opt.scale)
            # [(self.mutation,self.proliferation,self.size)]
            fpath = "{0}/{1}/{1}-circles_all.dat".format(self.opt.test_group_dir,
                                                         self.opt.param_set)
            mutation_v_proliferation_dat(self.subpop.tree_to_list("circles_all"),
                                         fpath,
                                         "MUTATION V PROLIFERATION RATES",
                                         self.opt.scale)


    def print_plots(self, when, end_time):
        """ Print all results to plots / file


        Print result to graphs using matplotlib
        If --R option parsed print raw output to file

        """

        filename = "{0}/{1}/{2}/{3}".format(self.opt.test_group_dir,
                                            self.opt.param_set,
                                            self.opt.run_number,
                                            when)

        if self.opt.R:
            from outputdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat
        else:
            from plotdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat, \
                    mutation_distribution, mutation_crash, \
                    make_subpop_life_mut, make_dual_hist, make_dual_box

        # Proliferation Histogram #

        if self.tumoursize > 0:
            end_proliferation = self.subpop.tree_to_list("proliferation_size")
            end_mutation = self.subpop.tree_to_list("mutation_rate")

            #PROLIFERATION HISTOGRAM
            # [(self.proliferation-self.prolif_adj,self.size)]

            end_proliferation.sort()
            self.mid_proliferation.sort()
            make_dual_hist(end_proliferation, self.mid_proliferation,
                           filename + "prolif_hist", "Proliferation Rates")

            end_mutation.sort()
            self.mid_mutation.sort()
            make_dual_hist(end_mutation, self.mid_mutation,
                           filename + "mutation_hist", "Mutation Rates")

            make_dual_box(end_mutation, self.mid_mutation,
                          filename + "mutation_box", "Mutation Rate Box")
