"""
Module containing the Treatment class.

This module contains the main
treatment class, and any treatment-
or chemo-related functions.

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
import tree_to_xml
import plotdata
import dropdata

class Treatment(object):
    """
    Class for modelling tumour treatment.

    This class determines when to introduce
    selective pressure; tracks the amount of
    'drug' in the tumour; and determines if
    and when to repeat a dose.
    """
    def __init__(self, opt):
        """Create a new Treatment object"""
        self.select_time = opt.select_time
        self.select_pressure = opt.select_pressure
        self.mutagenic_pressure = opt.mutagenic_pressure
        self.M = opt.M
        self.max_size_lim = opt.max_size_lim
        self.is_introduced = False

    def introduce(self, popn, t):
        """Introduce treatment into simulation."""
        self.is_introduced = True
        # note the time that treatment was introduced
        self.select_time = t
        # update population in various ways
        popn.select_pressure = self.select_pressure
        popn.mutagenic_pressure = self.mutagenic_pressure
        popn.subpop.set_precrash_size()
        popn.mid_proliferation = popn.subpop.tree_to_list("proliferation_size")
        popn.mid_mutation = popn.subpop.tree_to_list("mutation_rate")
        if not popn.opt.NP:
            plotdata.print_results(popn, "mid", t)
        tree_to_xml.tree_parse(popn.subpop, popn.tumoursize,
                               t, "mid0")
        if popn.opt.init_diversity:
            dropdata.drop(popn.subpop, popn.tumoursize,
                          t, "mid0")
        popn.selective_pressure_applied = True

    def update(self, popn, t):
        """Update treatment status."""
        if not self.is_introduced:
            if t == self.select_time:
                self.introduce(popn, t)
            elif self.M:
                # auto introduction of treatment
                if popn.exceeds_size_limit(self.max_size_lim):
                    self.introduce(popn, t)
