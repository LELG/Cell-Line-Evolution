"""
Tumour-level logic and data.

Authors
-------
Andrew Bakshi : andrew.bakshi@gmail.com
Yoshua Wakeham : yoshua.wakeham@petermac.org
"""

from __future__ import print_function
from analytics import Analytics
from subpopulation import Subpopulation


class Population(object):
    """
    Tumour population.

    This object tracks tumour-level data
    such as tumour size, clone count, and
    average mutation and proliferation rates.

    Its most significant attribute is subpop,
    a Subpopulation object, which functions as
    the root of the 'tree' of clones.
    """
    def __init__(self, opt, from_file=False):
        self.opt = opt
        self.max_size_lim = opt.max_size_lim
        self.prolif_lim = self.opt.prolif_lim

        if not from_file:
            self.tumoursize = opt.init_size
            self.clonecount = 1
            self.all_mutations = {'b': [], 'n': [], 'd': [], 'r': []}
            self.subpop = Subpopulation(opt=opt,
                                        prolif=opt.pro, mut_rate=opt.mut,
                                        depth=0, t_curr=0,
                                        col='n', prev_time=0)
            if opt.init_diversity:
                self.subpop.size = 0
                self.subpop.new_subpop_from_file(self.opt, self.opt.sub_file)
            else:
                self.subpop.size = opt.init_size
            self.avg_pro_rate = self.opt.pro
            self.avg_mut_rate = self.opt.mut
        else:
            # these attributes will be initialised from file
            self.tumoursize = self.clonecount = None
            self.all_mutations = self.subpop = None
            self.avg_pro_rate = self.avg_mut_rate = None

        self.analytics_base = Analytics()
        self.select_pressure = 0.0
        self.mutagenic_pressure = 0.0
        self.selective_pressure_applied = False
        # these lists will be populated at crash time
        self.mid_proliferation = []
        self.mid_mutation = []

    @classmethod
    def init_from_file(cls, opt, extra_params, root_clone, all_muts):
        """
        Initialise a population from a stored snapshot.

        Inputs
        ------
        opt: global parameter set, just like one created for
             a normal simulation
        extra_params: dictionary containing values for various
                      population attributes. See snapshot.py
                      for details.

        Returns
        -------
        A new Population object.
        """
        new_popn = cls(opt, from_file=True)
        new_popn.subpop = root_clone
        new_popn.all_mutations = all_muts

        # Update attributes from file.
        # Only update params that the population already has,
        # to prevent bad input from pickled file
        for param in extra_params:
            if hasattr(new_popn, param):
                setattr(new_popn, param, extra_params[param])
            else:
                raise Exception("attempting to set invalid attribute for population: {}".format(param))

    def update(self, treatmt, t_curr):
        """
        Update tumour and clones for a single time step.

        Check quantity of selective pressure present,
        recalculate ratio of population size to max size,
        update subclones, and recalculate average
        mutation and proliferation rates.

        Args
        ----
        treatmt : a treatment object
        t_curr : current time step

        Returns
        -------
        None
        """
        if treatmt.is_introduced:
            # update amount of selective and mutagenic pressure
            self.select_pressure = treatmt.curr_select_pressure
            self.mutagenic_pressure = treatmt.curr_mut_pressure

        # determine proliferation adjustment due to
        # population size
        ratio = self.tumoursize / float(self.max_size_lim)
        prolif_adj = ratio * self.prolif_lim

        if self.opt.prune_clones:
            # delete all dead, childless clones
            self.subpop.prune_dead_end_clones()

        # update subpopulations, getting back
        # tumour size, clone count, and aggregate
        # mutation and proliferation rates
        subpop_results = self.subpop.update(self.opt,
                                            self.select_pressure,
                                            self.mutagenic_pressure,
                                            t_curr, prolif_adj,
                                            self.all_mutations)
        self.tumoursize, self.clonecount, agg_mut, agg_pro = subpop_results

        if not self.is_dead():
            self.avg_mut_rate = agg_mut / float(self.tumoursize)
            self.avg_pro_rate = agg_pro / float(self.tumoursize)

    def is_dead(self):
        """Determine if this population has died out."""
        return self.tumoursize <= 0

    def exceeds_size_limit(self, max_size_lim, tolerance=0.0):
        """Determine whether this population has exceeded a certain size."""
        size_limit = max_size_lim * (1 + tolerance)
        return self.tumoursize > size_limit

    def record_treatment_introduction(self, treatmt):
        """Set and note certain variables at time of crash."""
        self.select_pressure = treatmt.curr_select_pressure
        self.mutagenic_pressure = treatmt.curr_mut_pressure
        self.selective_pressure_applied = True
        self.subpop.set_precrash_size()
        self.mid_proliferation = self.subpop.get_clone_attrs_as_list(["prolif_rate", "size"])
        self.mid_mutation = self.subpop.get_clone_attrs_as_list(["mut_rate", "size"])
