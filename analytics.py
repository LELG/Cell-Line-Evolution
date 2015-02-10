"""
Module defining classes and functions
related to simulation results and analysis.

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


class Analytics(object):
    """
    Simple class for recording analytics on a population.

    population
        population size at each point in time
    subpopulation
        number of subpopulations (clones) at each time
    time
        list of time 1... 1000 for plot purposes
    mutation
        Store average mutation rate at each point in time
    proliferation
        Store average proliferation at each point in time
    subpop_mutation
        Store mutation rate for each subpopulation
    subpop_proliferation
        Store Proliferation rate for each subpopulation
    """
    def __init__(self):
        self.population = []
        self.subpopulation = []
        self.time = []
        self.mutation = []
        self.proliferation = []
        #self.subpop_mutation = []
        #self.subpop_proliferation = []

    def update(self, popn, treatmt, t):
        """
        Update analytics for a single time step.
        """
        self.population.append(popn.tumoursize)
        self.subpopulation.append(popn.clonecount)
        self.time.append(t)

        # append effective proliferation
        if treatmt.is_introduced:
            eff_prolif = popn.avg_pro_rate - treatmt.select_pressure
            self.proliferation.append(eff_prolif)
        else:
            self.proliferation.append(popn.avg_pro_rate)

        self.mutation.append(popn.avg_mut_rate)


def went_through_crash(sim, popn):
    """Determine whether a population survived the crash."""
    CRASH_BUFFER = 25 #check just past crash time
    post_crash_time = sim.select_time + CRASH_BUFFER
    return len(popn.analytics_base.time) > post_crash_time


def completion_status(sim, popn):
    """
    Return recovery, full or partial,  percent

    eg. 'N'
    eg. 'Y'
    """
    recover = 'N'
    recover_type = 'NONE'

    if went_through_crash(sim, popn):
        if popn.tumoursize > sim.max_size_lim * 0.5:
            recover = 'Y'
            recover_type = 'PART'
            if popn.tumoursize > sim.max_size_lim * 0.75:
                # recovered population is > 75% of original size
                recover_type = 'FULL'
    else: #didnt crash
        if popn.tumoursize > sim.max_size_lim * 0.75:
            recover_type = 'FULLNC'
            recover = 'Y'

    recover_percent = popn.tumoursize / float(sim.max_size_lim)
    return recover, recover_type, recover_percent


def postcrash_minmax(sim, popn):
    """Return data on max and min post-crash population size."""
    post_crash_pop = popn.analytics_base.population[sim.select_time:]
    min_val = min(post_crash_pop) #VALIDATION - can return empty
    min_val_index = post_crash_pop[post_crash_pop.index(min_val)]
    #VALIDATION - can return empty
    #maximum value after lowest point
    max_val = 0
    max_time = 0
    if post_crash_pop[min_val_index:]:
        max_val = max(post_crash_pop[min_val_index:])
        max_time = post_crash_pop.index(max_val) + sim.select_time
    #add back time that we cut out when filtering post crash
    min_time = post_crash_pop.index(min_val) + sim.select_time
    #max time after low point
    return min_val, min_time, max_val, max_time


def precrash_minmax(sim, popn):
    pre_crash_pop = popn.analytics_base.population[:sim.select_time]
    min_val = min(pre_crash_pop)
    max_val = max(pre_crash_pop)
    min_time = pre_crash_pop.index(min_val)
    max_time = pre_crash_pop.index(max_val)
    return min_val, min_time, max_val, max_time
