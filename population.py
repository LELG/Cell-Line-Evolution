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

    Its most significant attribute is its
    Subpopulation, which functions as the root
    of the 'tree' of subclones.
    """
    def __init__(self, opt):
        self.opt = opt
        self.tumoursize = opt.init_size
        self.clonecount = 1
        self.max_size_lim = opt.max_size_lim
        self.prolif_lim = self.opt.prolif_lim
        self.subpop = Subpopulation(opt=vars(opt),
                                    p=opt.pro, m=opt.mut,
                                    depth=0, time=0,
                                    mut_type='n', col='n',
                                    pmp=opt.prob_mut_pos,
                                    pmn=opt.prob_mut_neg,
                                    pim=opt.prob_inc_mut,
                                    pdm=opt.prob_dec_mut,
                                    msc=opt.mscale, prev_time=0)
        if opt.init_diversity:
            self.subpop.size = 0 #don't use init dummy popn if reading from file
            self.subpop.newsubpop_from_file(self.opt.sub_file)
        else:
            self.subpop.size = opt.init_size
        self.analytics_base = Analytics()
        self.select_pressure = 0.0
        self.mutagenic_pressure = 0.0
        self.selective_pressure_applied = False
        self.avg_pro_rate = self.opt.pro
        self.avg_mut_rate = self.opt.mut
        # these lists will be populated at crash time
        self.mid_proliferation = []
        self.mid_mutation = []

    def update(self, treatmt, t):
        """
        Update tumour and clones for a single time step.

        Check quantity of selective pressure present,
        recalculate ratio of population size to max size,
        update subclones, and recalculate average
        mutation and proliferation rates.

        Args
        ----
        treatmt : a treatment object
        t : current time step

        Returns
        -------
        None
        """
        if treatmt.is_introduced:
            # update amount of selective and mutagenic pressure
            self.select_pressure = treatmt.select_pressure
            self.mutagenic_pressure = treatmt.mutagenic_pressure

        # determine proliferation adjustment due to
        # population size
        ratio = self.tumoursize / float(self.max_size_lim)
        prolif_adj = ratio * self.prolif_lim

        # update subpopulations, getting back
        # tumour size, clone count, and aggregate
        # mutation and proliferation rates
        subpop_results = self.subpop.cycle(self.tumoursize,
                                           self.select_pressure,
                                           self.mutagenic_pressure,
                                           t, prolif_adj)
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
