"""
High-level simulation routines and logic.

The Simulator object controls all high-level
aspects of the simulation, and co-ordinates
interactions between the Population (tumour)
and treatment objects.

Authors
-------
Yoshua Wakeham : yoshwakeham@gmail.com

Date created
------------
10 Feb 2015

Change log
----------
11 Feb 2015 - Various small changes
            - Add documentation
"""
import time
import os
import csv
import popln.population as population
import popln.treatment as treatment
import popln.analytics as analytics
from popln.utilities import secs_to_hms
import popln.tree_to_xml as tree_to_xml
import popln.plotdata as plotdata
import popln.dropdata as dropdata

END_POP_TOO_LARGE = "Population exceeded size limit."
END_POP_DIED_OUT = "Population died out."
END_MAX_CYCLES = "Simulation reached maximum cycle limit."

class Simulator(object):
    """
    Simulation controller.

    Control all high-level behaviour of the simulation,
    including running the simulation and co-ordinating
    interactions between the tumour and treatment objects.

    Attributes
    ----------
    opt : namespace containing all command-line
        parameters used to run this simulation.
    is_running : bool, record whether the sim
        is currently running.
    runtime : simulation runtime.
    popn_recovered : whether the population has
        survived a treament-induced crash.
    total_cycles : total number of time steps
        elapsed in simulation.
    popn : the tumour
    treatmt : the treatment
    """
    def __init__(self, opt):
        """
        Initialise Simulator from command line parameters.

        Note
        ----
        See popln.main.parse_cmd_line_args() for a
        description of the command line arguments.
        """
        self.opt = opt

        # calculate and record proliferation limit
        # this avoids the need to take it from command line
        self.opt.prolif_lim = opt.pro - opt.die

        # TODO determine which of these params
        # really need to be copied out of opt.
        self.test_group = opt.test_group
        self.param_set = opt.param_set
        self.run_number = opt.run_number

        self.test_group_dir = opt.test_group_dir
        self.param_set_dir = opt.param_set_dir
        self.run_dir = opt.run_dir

        self.max_cycles = opt.max_cycles
        self.max_size_lim = opt.max_size_lim

        self.init_diversity = opt.init_diversity
        self.sub_file = opt.sub_file

        # TODO shift these parameters to Treatment class?
        self.select_time = opt.select_time
        self.select_pressure = opt.select_pressure
        self.mutagenic_pressure = opt.mutagenic_pressure

        self.scale = opt.scale
        self.mscale = opt.mscale

        # TODO rename these options to something more meaningful
        self.R = opt.R
        self.M = opt.M
        self.Z = opt.Z
        self.NP = opt.NP

        self.is_running = False
        self.treatment_introduced = False
        self.runtime = None
        self.popn_recovered = False
        self.total_cycles = self.max_cycles

        # create the Population
        self.popn = population.Population(self.opt)
        # create the Treatment
        self.treatmt = treatment.Treatment(self.opt)


    def run(self):
        """
        Run simulation.

        Start the simulation, iterate over
        time steps, ending the simulation early
        if any of the end conditions are met. Also
        record the simulation runtime.

        Args
        ----
        None

        Returns
        -------
        None. Concludes with a call to Simulator.finish().
        """
        self.is_running = True

        # begin timing simulation
        start_time = time.time()

        end_condition = END_MAX_CYCLES
        for t in xrange(self.max_cycles):
            self.update(t)

            # test for end conditions
            if self.popn.exceeds_size_limit(self.max_size_lim, tolerance=0.05):
                if self.treatmt.is_introduced:
                    end_condition = END_POP_TOO_LARGE
                    self.total_cycles = t
                    self.popn_recovered = True
                    break
            if self.popn.is_dead():
                end_condition = END_POP_DIED_OUT
                self.total_cycles = t
                break

        # finish timing
        end_time = time.time()
        self.runtime = end_time - start_time

        # end simulation
        self.finish(end_condition)

    def update(self, t):
        """
        Simulate a single time step.

        Update the treatment (for instance,
        recalculating amount of drug remaining),
        population (reproduce/mutate/die), and
        analytics. Print a status update at
        specified time increments.

        Args
        ----
        t : Current time step.

        Returns
        -------
        None.
        """
        self.treatmt.update(self.popn, t)
        self.popn.update(self.treatmt, t)
        self.popn.analytics_base.update(self.popn, self.treatmt, t)

        # print status message
        if t % 1000 == 0:
            self.print_status_update(t)

    def finish(self, end_condition):
        print("SIMULATION ENDED: {}".format(end_condition))
        self.write_summary(self.popn, self.treatmt,
                           self.total_cycles, self.runtime,
                           self.popn_recovered)
        # test analytics dump function
        fpath = "{0}/alldata.csv".format(self.run_dir)
        self.popn.analytics_base.write_to_file(fpath)
        if not self.NP:
            plotdata.print_results(self.popn, "end", self.total_cycles)
            plotdata.print_plots(self.popn, "new")
        fname = ""
        tree_to_xml.tree_parse(self.popn.subpop, self.popn.tumoursize, self.total_cycles, fname)
        if self.init_diversity:
            print("Printing drop data")
            dropdata.drop(self.popn.subpop, self.popn.tumoursize, self.total_cycles, "end")

        self.is_running = False


    def print_status_update(self, t):
        os.system('clear')
        status_msg = """
Tumour Evolution Simulation
---------------------------

      ... running ...

Cycle:       {0} of {1}
Tumour size: {2}
"""
        print(status_msg.format(t, self.max_cycles, self.popn.tumoursize))

    def print_info(self):
        """Print information about this simulation."""
        print("Simulation Parameter Set")
        print("------------------------")
        for attr, val in self.opt:
            print("{}: {}".format(attr, val))

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

        tg_summary_fpath = "{0}/{1}_results.csv".format(self.test_group_dir,
                                                        self.test_group)
        ps_summary_fpath = "{0}/{1}_{2}_results.csv".format(self.param_set_dir,
                                                            self.test_group,
                                                            self.param_set)
        tg_results_file = open(tg_summary_fpath, 'a')
        ps_results_file = open(ps_summary_fpath, 'a')
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
