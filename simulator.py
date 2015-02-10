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
import population
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
    and treatment objects.
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
            if self.population_too_big(self.popn) and self.treatment_introduced:
                end_condition = POP_TOO_LARGE
                self.total_cycles = t
                self.popn_recovered = True
                break
            if self.population_died_out(self.popn):
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
        self.popn.cycle(t)

        if not self.treatment_introduced:
            if t == self.select_time:
                self.popn.subpop.set_precrash_size()
                self.popn.selective_pressure()
                if not self.NP:
                    self.popn.print_results("mid", t)
                tree_to_xml.tree_parse(self.popn.subpop, self.popn.tumoursize,
                                       t, "mid0")
                if self.init_diversity:
                    dropdata.drop(self.popn.subpop, self.popn.tumoursize,
                                  t, "mid0")
                self.treatment_introduced = True
                self.popn.selective_pressure_applied = True
            elif self.M:
                # auto introduction of treatment
                if self.population_too_big(self.popn):
                    self.popn.subpop.set_precrash_size()
                    self.popn.selective_pressure()
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
                    self.treatment_introduced = True
                    self.popn.selective_pressure_applied = True
        # print status message
        if t % 1000 == 0:
            self.print_status_update(t)

    def finish(self, end_condition):
        print("Simulation complete: {}".format(end_condition))
        self.popn.write_population_summary(self.total_cycles, self.runtime, self.popn_recovered)
        if not self.NP:
            self.popn.print_results("end", self.total_cycles)
            self.popn.print_plots("new", self.total_cycles)
        fname = ""
        tree_to_xml.tree_parse(self.popn.subpop, self.popn.tumoursize, self.total_cycles, fname)
        if self.init_diversity:
            print("Printing drop data")
            dropdata.drop(self.popn.subpop, self.popn.tumoursize, self.total_cycles, "end")

    def population_too_big(self, popn):
        SIZE_TOLERANCE = 0.05
        tumour_size_threshold = self.max_size_lim * (1 + SIZE_TOLERANCE)
        return popn.tumoursize > tumour_size_threshold

    def population_died_out(self, popn):
        return popn.tumoursize <= 0

    def print_status_update(self, t):
        os.system('clear')
        print("Cycle {0} of {1}: Tumour size: {2}".format(t, self.max_cycles,
                                                          self.popn.tumoursize))
