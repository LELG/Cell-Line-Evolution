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
        # adopt relevant command line parameters
        self.treatment_type = opt.treatment_type
        self.decay_type = opt.decay_type
        self.decay_rate = opt.decay_rate
        self.decay_func = None
        self.select_time = opt.select_time
        self.init_select_pressure = opt.select_pressure
        self.init_mut_pressure = opt.mutagenic_pressure

        self.auto_treatment = opt.auto_treatment
        self.max_size_lim = opt.max_size_lim

        # assume that all simulations will start
        # with no selective pressure introduced
        self.curr_select_pressure = 0.0
        self.curr_mut_pressure = 0.0

        self.is_introduced = False


    def introduce(self, popn, t_curr):
        """
        Introduce treatment into simulation.

        Let the tumour know that the treatment
        has been introduced, and produce some plots.
        """
        self.is_introduced = True
        self.select_time = t_curr
        self.curr_select_pressure = self.init_select_pressure
        self.curr_mut_pressure = self.init_mut_pressure
        # set decay function
        self.decay_func = get_decay_func(self.decay_type,
                                         self.decay_rate,
                                         self.init_select_pressure,
                                         t_curr)
        # let population know that treatment has
        # been introduced
        popn.record_treatment_introduction(self)
        # TODO move this plotting elsewhere
        if not popn.opt.no_plots:
            plotdata.print_results(popn, "mid", t_curr)
        tree_to_xml.tree_parse(popn.subpop, popn.tumoursize,
                               t_curr, "mid0")
        if popn.opt.init_diversity:
            dropdata.drop(popn.subpop, popn.tumoursize,
                          t_curr, "mid0")

    def update(self, popn, t_curr):
        """
        Update treatment status.

        If treatment has not yet been introduced,
        check whether (any of) the conditions for
        introduction are now met.
        """
        if self.is_introduced:
            self.decay(t_curr)
            if self.reintroduction_conditions_met(popn, t_curr):
                self.reintroduce(popn, t_curr)
        else:
            # not yet introduced
            if t_curr == self.select_time:
                self.introduce(popn, t_curr)
            elif self.auto_treatment:
                if popn.exceeds_size_limit(self.max_size_lim):
                    self.introduce(popn, t_curr)

    def decay(self, t_curr):
        """Model selective pressure decay for this time step."""
        if self.sel_pressure_remaining():
            self.curr_select_pressure = self.decay_func(t_curr)

    def sel_pressure_remaining(self):
        """Determine if any pressure remains in the tumour."""
        return self.curr_select_pressure > 0.0

    def reintroduction_conditions_met(self, popn, t_curr):
        pass

    def reintroduce(self, popn, t_curr):
        pass


class MetronomicTreatment(Treatment):
    """Model metronomic treatment, that is, small doses at frequent, regular intervals."""
    def __init__(self, opt):
        super(MetronomicTreatment, self).__init__(opt)
        self.cycles_per_dose = 100
        # scale down the treatment to model metronomic treatment
        # self.init_select_pressure *= 0.5
        # self.init_mut_pressure *= 0.5

    def reintroduction_conditions_met(self, popn, t_curr):
        t_delta = t_curr - self.select_time
        return t_delta >= self.cycles_per_dose

    def reintroduce(self, popn, t_curr):
        self.curr_select_pressure += self.init_select_pressure
        self.curr_mut_pressure += self.init_mut_pressure
        self.select_time = t_curr
        # reset decay function
        self.decay_func = get_decay_func(self.decay_type,
                                         self.decay_rate,
                                         self.curr_select_pressure,
                                         t_curr)


def get_decay_func(decay_type, decay_rate, init_qty, t_init):
    """Return a partially-applied decay function that will take a single time parameter."""
    if decay_type == 'constant':
        return partial(constant_decay, init_qty)
    elif decay_type == 'linear':
        return partial(linear_decay, init_qty, decay_rate, t_init)
    elif decay_type == 'exp':
        return partial(exp_decay, init_qty, decay_rate, t_init)


def constant_decay(init_qty, t_curr):
    """Model a constant (non-decaying) quantity."""
    return init_qty


def linear_decay(init_qty, decay_rate, t_init, t_curr):
    """Model a linearly decaying quantity."""
    t_delta = t_curr - t_init
    return init_qty - decay_rate * t_delta


def exp_decay(init_qty, decay_rate, t_init, t_curr):
    """Model an exponentially decaying quantity."""
    t_delta = t_curr - t_init
    return init_qty * math.exp(-1 * decay_rate * t_delta)
