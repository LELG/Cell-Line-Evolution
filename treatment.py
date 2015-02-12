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
import math
from functools import partial

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

        #self.decay_type = 'constant'
        #self.decay_type = 'linear'
        #self.decay_rate = 0.00002
        self.decay_type = 'exp'
        self.decay_rate = 0.001
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
        # set decay function
        self.decay_func = get_decay_func(self.decay_type,
                                         self.decay_rate,
                                         self.init_select_pressure,
                                         t)
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


def get_decay_func(decay_type, decay_rate, init_qty, t_init):
    if decay_type == 'constant':
        return partial(constant_decay, init_qty)
    elif decay_type == 'linear':
        return partial(linear_decay, init_qty, decay_rate, t_init)
    elif decay_type == 'exp':
        return partial(exp_decay, init_qty, decay_rate, t_init)


def constant_decay(init_qty, t_curr):
    return init_qty


def linear_decay(init_qty, decay_rate, t_init, t_curr):
    t_delta = t_curr - t_init
    return init_qty - decay_rate * t_delta


def exp_decay(init_qty, decay_rate, t_init, t_curr):
    t_delta = t_curr - t_init
    return init_qty * math.exp(-1 * decay_rate * t_delta)
