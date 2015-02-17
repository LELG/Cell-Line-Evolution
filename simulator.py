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
from textwrap import dedent
import population
import treatment
import analytics
from utilities import secs_to_hms
import tree_to_xml
import plotdata
import dropdata

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

        self.r_output = opt.r_output
        self.auto_treatment = opt.auto_treatment
        self.prune_clones = opt.prune_clones
        self.no_plots = opt.no_plots

        self.is_running = False
        self.treatment_introduced = False
        self.runtime = None
        self.total_cycles = self.max_cycles

        # create the Population
        self.popn = population.Population(self.opt)
        # create the Treatment
        if opt.treatment_type == 'single_dose':
            self.treatmt = treatment.Treatment(self.opt)
        elif opt.treatment_type == 'metronomic':
            self.treatmt = treatment.MetronomicTreatment(self.opt)
        elif opt.treatment_type == 'adaptive':
            self.treatmt = treatment.AdaptiveTreatment(self.opt)
        else:
            raise ValueError("Bad value for treatment type parameter")


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
                    break
            if self.popn.is_dead():
                end_condition = END_POP_DIED_OUT
                self.total_cycles = t
                break

        # finish timing
        end_time = time.time()
        self.runtime = end_time - start_time

        self.is_running = False

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
        """
        Complete the simulation.

        Write summary files and print plots
        to complete this sim run.

        Args
        ----
        end_condition : a human-readable message
            describing why this simulation ended

        Returns
        -------
        None.
        """
        print("SIMULATION ENDED: {}".format(end_condition))
        # write to summary file
        self.write_summary(self.popn, self.treatmt,
                           self.total_cycles, self.runtime)
        # dump all run data to CSV file
        data_dump_fpath = "{0}/{1}_{2}_{3}_alldata.csv".format(self.run_dir,
                                                               self.test_group,
                                                               self.param_set,
                                                               self.run_number)
        self.popn.analytics_base.write_to_file(data_dump_fpath)
        # make plots
        if not self.no_plots:
            plotdata.print_results(self.popn, "end", self.total_cycles)
            plotdata.print_plots(self.popn, "new")
        fname = ""
        # write phylogenetic tree to XML file
        tree_to_xml.tree_parse(self.popn.subpop, self.popn.tumoursize,
                               self.total_cycles, fname)
        # if heterogeneous initial pop, output drop data
        if self.init_diversity:
            print("Printing drop data")
            dropdata.drop(self.popn.subpop, self.popn.tumoursize,
                          self.total_cycles, "end")


    def print_status_update(self, t):
        """ Print current time step and tumour size for running sim."""
        os.system('clear')
        status_msg = dedent("""\
            Tumour Evolution Simulation
            ---------------------------

                  ... running ...

            Cycle:       {0} of {1}
            Tumour size: {2}
            """)
        print(status_msg.format(t, self.max_cycles, self.popn.tumoursize))

    def print_info(self):
        """Print simulation's initial parameter set."""
        hdr = dedent("""\
            Simulation Parameter Set
            ------------------------
            """)
        print(hdr)
        for attr, val in vars(self.opt).items():
            print("{}: {}".format(attr, val))
        print('\n')

    def write_summary(self, popn, treatmt, tot_cycles, elapsed_time):
        """
        Write simulation summary to file(s).

        Write summary data about this simulation run
        to two summary files, in CSV format:

        'testgroup_results.csv' (master summary file for this test group)
        'testgroup_paramset_results.csv' (summary for this param set only)

        Args
        ----
        popn : a tumour population
        treatmt : a treatment regime
        tot_cycles : total number of cycles for this simulation
        elapsed_time : simulation runtime in seconds

        Returns
        -------
        None
        """
        # Get pre-crash min and max population size, and
        # the times at which they occurred
        precrash_minmax_data = analytics.precrash_minmax(treatmt, popn)
        min_val, min_time, max_val, max_time = precrash_minmax_data

        cmin_val = cmin_time = cmax_val = cmax_time = 0

        went_through_crash = 'N'
        if analytics.went_through_crash(treatmt, popn):
            went_through_crash = 'Y'
            #Get post-crash min and max values
            postcrash_minmax_data = analytics.postcrash_minmax(treatmt, popn)
            cmin_val, cmin_time, cmax_val, cmax_time = postcrash_minmax_data
            #hasty fix for calculating max time
            # if cmax_time == 0 and recovered:
            #    cmax_time = tot_cycles

        # determine whether, when and how fully the population recovered
        recovery_status = analytics.completion_status(self, treatmt, popn)
        recovered, recover_type, recover_percent = recovery_status

        # open files (these should have been created already)
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
        summary_vals = (self.param_set, self.run_number, went_through_crash,
                        recovered, recover_type, recover_percent,
                        popn.opt.pro, popn.opt.die, popn.opt.mut,
                        treatmt.select_time, treatmt.init_select_pressure,
                        popn.opt.prob_mut_pos, popn.opt.prob_mut_neg,
                        popn.opt.prob_inc_mut, popn.opt.prob_dec_mut,
                        popn.analytics_base.tumoursize[-1],
                        popn.analytics_base.clonecount[-1],
                        popn.analytics_base.avg_mutation[-1],
                        popn.analytics_base.avg_proliferation[-1],
                        secs_to_hms(elapsed_time), tot_cycles,
                        min_val, min_time, max_val, max_time,
                        cmin_val, cmin_time, cmax_val, cmax_time)

        tg_writer.writerow(summary_vals)
        ps_writer.writerow(summary_vals)
        tg_results_file.close()
        ps_results_file.close()