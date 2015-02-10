"""
Class for controllong entire simulation.

The Simulator object controls all high-level
aspects of the simulation, and co-ordinates
interactions between the Population (tumour)
and treatment objects.

Author
------
Yoshua Wakeham : yoshwakeham@gmail.com

Date created
------------
10 Feb 2015

Notes
-----

Change log
----------
"""
import time
import os
import csv
import population
import treatment
import analytics
from utilities import secs_to_hms
import tree_to_xml
import dropdata

POP_TOO_LARGE = "Population exceeded size limit."
POP_DIED_OUT = "Population died out."
MAX_CYCLES = "Simulation reached maximum cycle limit."

class Simulator(object):
    """
    Master class representing entire simulation.

    The Simulator object controls all high-level
    aspects of the simulation, and co-ordinates
    interactions between the Population (tumour)
    and Treatment objects.
    """
    def __init__(self, parameters):
        """
        Create new Simulator from command line parameters.
        """
        self.test_group = parameters.test_group
        self.test_group_dir = parameters.test_group_dir
        self.param_set = parameters.param_set
        self.run_number = parameters.run_number

        self.max_cycles = parameters.max_cycles
        self.max_size_lim = parameters.max_size_lim

        self.init_diversity = parameters.init_diversity
        self.sub_file = parameters.sub_file

        # TODO shift these parameters to Treatment class
        self.select_time = parameters.select_time
        self.select_pressure = parameters.select_pressure
        self.mutagenic_pressure = parameters.mutagenic_pressure

        self.scale = parameters.scale
        self.mscale = parameters.mscale

        self.R = parameters.R
        self.M = parameters.M
        self.Z = parameters.Z
        self.NP = parameters.NP

        self.is_running = False
        self.treatment_introduced = False
        self.runtime = None
        self.popn_recovered = False
        self.total_cycles = self.max_cycles

        # create a Population object
        self.popn = population.Population(parameters)
        # create a Treatment object
        self.treatmt = treatment.Treatment(parameters)


    def run(self):
        """
        Run the simulation.
        """
        self.is_running = True

        # begin timing simulation
        start_time = time.time()

        end_condition = MAX_CYCLES
        for t in xrange(self.max_cycles):
            self.update(t)
            if self.popn.exceeds_size_limit(self.max_size_lim, tolerance=0.05):
                if self.treatmt.is_introduced:
                    end_condition = POP_TOO_LARGE
                    self.total_cycles = t
                    self.popn_recovered = True
                    break
            if self.popn.is_dead():
                end_condition = POP_DIED_OUT
                self.total_cycles = t
                break

        # finish timing
        end_time = time.time()
        self.runtime = end_time - start_time

        self.finish(end_condition)

    def update(self, t):
        """
        Simulate a single time step.
        """
        self.treatmt.update(self.popn, t)
        self.popn.update(self.treatmt, t)
        self.popn.analytics_base.update(self, self.popn, t)

        """
        if not self.treatmt.is_introduced:
            if t == self.select_time:
                self.popn.subpop.set_precrash_size()
                self.treatmt.introduce(self.popn, t)
                if not self.NP:
                    self.popn.print_results("mid", t)
                tree_to_xml.tree_parse(self.popn.subpop, self.popn.tumoursize,
                                       t, "mid0")
                if self.init_diversity:
                    dropdata.drop(self.popn.subpop, self.popn.tumoursize,
                                  t, "mid0")
                # TODO deprecate these flags; should be unnecessary
                self.treatment_introduced = True
                self.popn.selective_pressure_applied = True
            elif self.treatmt.M:
                # auto introduction of treatment
                if population_too_big(self.popn, self.max_size_lim):
                    self.popn.subpop.set_precrash_size()
                    self.treatmt.introduce(self.popn, t)
                    if not self.NP:
                        self.popn.print_results("mid", t)
                    tree_to_xml.tree_parse(self.popn.subpop,
                                           self.popn.tumoursize,
                                           t, "mid")
                    if self.init_diversity:
                        dropdata.drop(self.popn.subpop, self.popn.tumoursize,
                                      t, "mid")
                    #update time when sel press introduced
                    self.select_time = t

                    # TODO deprecate these, as above
                    self.treatment_introduced = True
                    self.popn.selective_pressure_applied = True
        """
        # print status message
        if t % 1000 == 0:
            self.print_status_update(t)

    def finish(self, end_condition):
        print("SIMULATION ENDED: {}".format(end_condition))
        self.write_summary(self.popn, self.treatmt,
                           self.total_cycles, self.runtime,
                           self.popn_recovered)
        if not self.NP:
            self.popn.print_results("end", self.total_cycles)
            self.popn.print_plots("new", self.total_cycles)
        fname = ""
        tree_to_xml.tree_parse(self.popn.subpop, self.popn.tumoursize, self.total_cycles, fname)
        if self.init_diversity:
            print("Printing drop data")
            dropdata.drop(self.popn.subpop, self.popn.tumoursize, self.total_cycles, "end")

        self.is_running = False


    def print_status_update(self, t):
        #os.system('clear')
        status_msg = """
Tumour Evolution Simulation
---------------------------

      ... running ...

Cycle:       {0} of {1}
Tumour size: {2}
"""
        print(status_msg.format(t, self.max_cycles, self.popn.tumoursize))

    def write_summary(self, popn, treatmt, num_cycles, elapsed_time, recovered):
        """
        Write simulation summary to file(s).

        Write summary data about this simulation run
        to two summary files, in CSV format:

        'testgroup_results.csv' (master summary file for this test group)
        'testgroup_paramset_results.csv' (summary for this param set only)
        """

        # localtime = time.asctime( time.localtime(time.time()) )

        #Get min and max values pre crash
        min_val, min_time, max_val, max_time = analytics.precrash_minmax(self, popn)
        cmin_val = cmin_time = cmax_val = cmax_time = 0

        went_through_crash = 'N'
        if analytics.went_through_crash(self, popn):
            went_through_crash = 'Y'
            #Get min and max values post crash
            #if survived past the crash + buffer
            cmin_val, cmin_time, cmax_val, cmax_time = analytics.postcrash_minmax(self, popn)
            #hasty fix for calculting max time
            if cmax_time == 0 and recovered:
                cmax_time = num_cycles

        recover, recover_type, recover_percent = analytics.completion_status(self, popn)

        tg_res_fpath = "{0}/{1}_results.csv".format(self.test_group_dir,
                                                    self.test_group)
        ps_res_fpath = "{0}/{1}/{2}_{1}_results.csv".format(self.test_group_dir,
                                                            self.param_set,
                                                            self.test_group)
        tg_results_file = open(tg_res_fpath, 'a')
        ps_results_file = open(ps_res_fpath, 'a')
        tg_writer = csv.writer(tg_results_file)
        ps_writer = csv.writer(ps_results_file)

        # assemble values to write
        summary_vals = (self.param_set, self.run_number,
                        went_through_crash,
                        recover, recover_type, recover_percent,
                        popn.opt.pro, popn.opt.die, popn.opt.mut,
                        treatmt.select_time, treatmt.select_pressure,
                        popn.opt.prob_mut_pos, popn.opt.prob_mut_neg,
                        popn.opt.prob_inc_mut, popn.opt.prob_dec_mut,
                        popn.analytics_base.population[-1],    # pop_size
                        popn.analytics_base.subpopulation[-1], # num_clones
                        popn.analytics_base.mutation[-1],      # avg_mut_rate
                        popn.analytics_base.proliferation[-1], # avg_pro_rate
                        secs_to_hms(elapsed_time),
                        num_cycles,
                        min_val, min_time,
                        max_val, max_time,
                        cmin_val, cmin_time,
                        cmax_val, cmax_time)

        tg_writer.writerow(summary_vals)
        ps_writer.writerow(summary_vals)
        tg_results_file.close()
        ps_results_file.close()

"""
    def selective_pressure(self, popn):
        popn.select_pressure = self.select_pressure
        popn.mutagenic_pressure = self.mutagenic_pressure
        popn.mid_proliferation = popn.subpop.tree_to_list("proliferation_size")
        popn.mid_mutation = popn.subpop.tree_to_list("mutation_rate")
"""
