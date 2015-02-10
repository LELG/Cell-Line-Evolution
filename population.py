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
from subpopulation import Subpopulation
import tree_to_xml
import time
import dropdata
from utilities import secs_to_hms


class Analytics(object):
    """
    Simple class for recording analytics on population.

    population
        population size at each point in time
    subpopulation
        number of subpopulations (clones) at each time
    time
        list of time 1... 1000 for plot purposes
    mutation
        Store average mutation rate at each point in time
    proliferation
        Store average proliferation at each point in time
    subpop_mutation
        Store mutation rate for each subpopulation
    subpop_proliferation
        Store Proliferation rate for each subpopulation
    """
    def __init__(self):
        self.population = []
        self.subpopulation = []
        self.time = []
        self.mutation = []
        self.proliferation = []
        self.subpop_mutation = []
        self.subpop_proliferation = []

    def update(self, population, t, avg_mut_rate, avg_pro_rate):
        self.population.append(population.tumoursize)
        self.subpopulation.append(population.clonecount)
        self.time.append(t)

        #EFFECTIVE PROLIFERATION
        if self.time[-1] > population.opt.select_time:
            self.proliferation.append(avg_pro_rate - population.opt.select_pressure)
        else:
            self.proliferation.append(avg_pro_rate)

        self.mutation.append(avg_mut_rate)


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
        if opt.init_diversity:
            self.subpop.size = 0 #don't use init dummy popn if reading from file
            self.subpop.newsubpop_from_file(self.opt.sub_file)


    def write_population_summary(self, num_cycles, elapsed_time, recovered):
        """
        Write simulation summary to file(s).

        Write summary data about this simulation run
        to two summary files, in CSV format:

        'testgroup_results.csv' (master summary file for this test group)
        'testgroup_paramset_results.csv' (summary for this param set only)
        """

        # localtime = time.asctime( time.localtime(time.time()) )

        #Get min and max values pre crash
        min_val, min_time, max_val, max_time = self.precrash_minmax()
        cmin_val = cmin_time = cmax_val = cmax_time = 0

        if self.went_through_crash():
            #Get min and max values post crash
            #if survived past the crash + buffer
            cmin_val, cmin_time, cmax_val, cmax_time = self.postcrash_minmax()
            #hasty fix for calculting max time
            if cmax_time == 0 and recovered:
                cmax_time = num_cycles

        recover, recover_type, recover_percent = self.complete_status()

        went_through_crash = 'N'
        #size_from_precrash = 0
        if self.went_through_crash():
            went_through_crash = 'Y'
            #pp = self.subpop.tree_to_list("two_side_size")
            #print("pp",pp)
            #for (pre, pos) in pp:
            #    if pre > 0:
            #        size_from_precrash += pos


        #TRACE self.opt.pro

        tg_res_fpath = "{0}/{1}_results.csv".format(self.opt.test_group_dir,
                                                    self.opt.test_group)
        ps_res_fpath = "{0}/{1}/{2}_{1}_results.csv".format(self.opt.test_group_dir,
                                                            self.opt.param_set,
                                                            self.opt.test_group)
        tg_results_file = open(tg_res_fpath, 'a')
        ps_results_file = open(ps_res_fpath, 'a')
        tg_writer = csv.writer(tg_results_file)
        ps_writer = csv.writer(ps_results_file)

        # assemble values to write
        summary_vals = (self.opt.param_set, self.opt.run_number,
                        went_through_crash,
                        recover, recover_type, recover_percent,
                        self.opt.pro, self.opt.die, self.opt.mut,
                        self.opt.select_time, self.opt.select_pressure,
                        self.subpop.prob_mut_pos, self.subpop.prob_mut_neg,
                        self.subpop.prob_inc_mut, self.subpop.prob_dec_mut,
                        self.analytics_base.population[-1],    # pop_size
                        self.analytics_base.subpopulation[-1], # num_clones
                        self.analytics_base.mutation[-1],      # avg_mut_rate
                        self.analytics_base.proliferation[-1], # avg_pro_rate
                        secs_to_hms(elapsed_time),
                        len(self.analytics_base.population),   # elapsed_cycles
                        min_val, min_time,
                        max_val, max_time,
                        cmin_val, cmin_time,
                        cmax_val, cmax_time)

        tg_writer.writerow(summary_vals)
        ps_writer.writerow(summary_vals)
        tg_results_file.close()
        ps_results_file.close()


    def went_through_crash(self):
        """Determine whether this population survived the crash."""
        crash_buffer = 25 #check just past crash time
        post_crash_time = self.opt.select_time + crash_buffer
        return len(self.analytics_base.population) > post_crash_time


    def complete_status(self):
        """
        Return recovery, full or partial,  percent

        eg. 'N'
        eg. 'Y'
        """
        recover = 'N'
        recover_type = 'NONE'

        if self.went_through_crash():
            if self.tumoursize > (self.max_size_lim/2):
                recover = 'Y'
                recover_type = 'PART'
                if self.tumoursize > self.max_size_lim * 0.75:
                    # recovered population is > 75% of original size
                    recover_type = 'FULL'
        else: #didnt crash
            if self.tumoursize > self.max_size_lim * 0.75:
                recover_type = 'FULLNC'
                recover = 'Y'

        recover_percent = self.tumoursize / float(self.max_size_lim)

        return recover, recover_type, recover_percent


    def postcrash_minmax(self):
        """Return data on max and min post-crash population size."""
        post_crash_pop = self.analytics_base.population[self.opt.select_time:]
        min_val = min(post_crash_pop) #VALIDATION - can return empty
        min_val_index = post_crash_pop[post_crash_pop.index(min_val)]
        #VALIDATION - can return empty
        #maximum value after lowest point
        max_val = 0
        max_time = 0
        if post_crash_pop[min_val_index:]:
            max_val = max(post_crash_pop[min_val_index:])
            max_time = post_crash_pop.index(max_val) + self.opt.select_time
        #add back time that we cut out when filtering post crash
        min_time = post_crash_pop.index(min_val) + self.opt.select_time
        #max time after low point
        return min_val, min_time, max_val, max_time

    def precrash_minmax(self):
        pre_crash_pop =\
                self.analytics_base.population[:self.opt.select_time]
        min_val = min(pre_crash_pop)
        max_val = max(pre_crash_pop)
        min_time = pre_crash_pop.index(min_val)
        max_time = pre_crash_pop.index(max_val)
        return min_val, min_time, max_val, max_time

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
        avg_mut_rate = 0
        avg_pro_rate = 0
        if self.tumoursize > 0:
            avg_mut_rate = agg_mut / float(self.tumoursize)
            avg_pro_rate = agg_pro / float(self.tumoursize)

        self.analytics_base.update(self, t, avg_mut_rate, avg_pro_rate)

    def selective_pressure(self):
        self.select_pressure = self.opt.select_pressure
        self.mutagenic_pressure = self.opt.mutagenic_pressure
        self.mid_proliferation = self.subpop.tree_to_list("proliferation_size")
        self.mid_mutation = self.subpop.tree_to_list("mutation_rate")


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
