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
    selective pressure present in the tumour
    at any given time step; and determines if
    and when to repeat a dose.

    Attributes
    ----------
    select_time : time step to introduce selective
        pressure, if other conditions have not already
        caused it to be introduced
    init_select_pressure : initial quantity of selective
        pressure to be introduced
    init_mut_pressure : initial quantity of mutagenic
        pressure to be introduced
    curr_select_pressure : amount of selective pressure
        remaining at this time step
    curr_mut_pressure : amount of mutagenic pressure remaining
        at this time step
    M : whether or not to introduce treatment automatically
        when population reaches a certain size
    max_size_lim : tumour size at which to introduce
        treatment automatically
    is_introduced : boolean, whether or not
        introduction has occurred
    decay_func : function to determine how selective pressure
        decays over time
    """
    def __init__(self, opt):
        """Create a new Treatment object"""
        self.select_time = opt.select_time
        self.init_select_pressure = opt.select_pressure
        self.init_mut_pressure = opt.mutagenic_pressure
        # this assumes that all simulations
        # will start with no drug introduced
        self.curr_select_pressure = 0.0
        self.curr_mut_pressure = 0.0
        self.M = opt.M
        self.max_size_lim = opt.max_size_lim
        self.is_introduced = False
        self.decay_func = None

    def introduce(self, popn, t):
        """
        Introduce treatment into simulation.

        Let the tumour know that the treatment
        has been introduced, and produce some plots.
        """
        self.is_introduced = True
        self.select_time = t
        self.curr_select_pressure = self.init_select_pressure
        self.curr_mut_pressure = self.init_mut_pressure
        self.decay_func = get_constant_decay(self.init_select_pressure)
        #self.decay_func = get_linear_decay(self.init_select_pressure,
        #                                   t, decay_rate=0.00002)
        # let population know that treatment has
        # been introduced
        popn.record_treatment_introduction(self)
        # TODO move this plotting elsewhere
        if not popn.opt.NP:
            plotdata.print_results(popn, "mid", t)
        tree_to_xml.tree_parse(popn.subpop, popn.tumoursize,
                               t, "mid0")
        if popn.opt.init_diversity:
            dropdata.drop(popn.subpop, popn.tumoursize,
                          t, "mid0")

    def update(self, popn, t):
        """
        Update treatment status.

        If treatment has not yet been introduced,
        check whether (any of) the conditions for
        introduction are now met.
        """
        if self.is_introduced:
            if self.sel_pressure_remaining():
                self.curr_select_pressure = self.decay_func(t)
        else:
            # not yet introduced
            if t == self.select_time:
                self.introduce(popn, t)
            elif self.M:
                # auto introduction of treatment
                if popn.exceeds_size_limit(self.max_size_lim):
                    self.introduce(popn, t)

    def sel_pressure_remaining(self):
        return self.curr_select_pressure > 0.0


def get_linear_decay(initial_quantity, t_init, decay_rate):
    """Return a linear decay function."""
    def linear_decay_func(t_curr):
        return initial_quantity - decay_rate * (t_curr - t_init)
    return linear_decay_func

def get_constant_decay(initial_quantity):
    """Return a constant decay (i.e. no decay) function."""
    return lambda t: initial_quantity
