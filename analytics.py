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
import csv

class Analytics(object):
    """
    Record analytics on a population.

    Attributes
    ----------

    population : population size at each time step
    subpopulation : number of clones at each time step
    time : time step (mainly tracked for plotting purposes)
    mutation : average mutation rate at each point in time
    proliferation : average proliferation rate at each point in time
    """
    def __init__(self):
        self.tumoursize = []
        self.clonecount = []
        self.time = []
        self.avg_mutation = []
        self.avg_proliferation = []
        self.select_pressure = []

    def update(self, popn, treatmt, t):
        """
        Update analytics for a single time step.

        Record the value of all tracked variables
        for this time step.
        """
        self.tumoursize.append(popn.tumoursize)
        self.clonecount.append(popn.clonecount)
        self.time.append(t)
        self.select_pressure.append(treatmt.curr_select_pressure)

        # append effective proliferation
        if treatmt.is_introduced:
            eff_avg_prolif = popn.avg_pro_rate - treatmt.curr_select_pressure
            self.avg_proliferation.append(eff_avg_prolif)
        else:
            self.avg_proliferation.append(popn.avg_pro_rate)

        self.avg_mutation.append(popn.avg_mut_rate)

    def write_to_file(self, filepath):
        """
        Write all analytics data to a CSV file.
        """
        data_file = open(filepath, 'w')
        writer = csv.writer(data_file)

        # get list of attributes
        data_dict = vars(self)
        # construct list of lists of the form
        # ['name', val0, val1, ... , valN]
        data_lists = [[varname] + data_dict[varname] for varname in data_dict]
        # transpose data using * magic
        data_transposed = zip(*data_lists)
        # write transposed data to file
        writer.writerows(data_transposed)
        data_file.close()


def went_through_crash(treatmt, popn):
    """Determine whether a population survived the crash."""
    CRASH_BUFFER = 25 #check just past crash time
    post_crash_time = treatmt.select_time + CRASH_BUFFER
    return len(popn.analytics_base.time) > post_crash_time


def completion_status(sim, treatmt, popn):
    """
    Return whether, and how fully, the population recovered from a crash.

    Returns
    -------
    recovered : whether or not the population survived the crash ('Y' or 'N')
    recover_type : whether the population recovered fully, or
        only partially
    recover_percent : ratio of final population size to size
        at time of crash
    """
    recovered = 'N'
    recover_type = 'NONE'

    if went_through_crash(treatmt, popn):
        if popn.tumoursize > sim.max_size_lim * 0.5:
            recovered = 'Y'
            recover_type = 'PART'
            if popn.tumoursize > sim.max_size_lim * 0.75:
                # recovered population is > 75% of original size
                recover_type = 'FULL'
    else: #didnt crash
        if popn.tumoursize > sim.max_size_lim * 0.75:
            recovered = 'Y'
            recover_type = 'FULLNC'

    recover_percent = popn.tumoursize / float(sim.max_size_lim)
    return recovered, recover_type, recover_percent


def postcrash_minmax(treatmt, popn):
    """
    Return data on max and min post-crash population size.

    Return the maximum and minimum post-crash population sizes,
    and the corresponding time steps.
    """
    post_crash_pop = popn.analytics_base.tumoursize[treatmt.select_time:]
    min_val = min(post_crash_pop) #VALIDATION - can return empty
    min_time = post_crash_pop.index(min_val) + treatmt.select_time
    #VALIDATION - can return empty
    #maximum value after lowest point
    max_val = 0
    max_time = 0
    recovering_pop = post_crash_pop[min_time:]
    if recovering_pop:
        max_val = max(recovering_pop)
        max_time = post_crash_pop.index(max_val) + treatmt.select_time
    return min_val, min_time, max_val, max_time


def precrash_minmax(treatmt, popn):
    """
    Return data on max and min pre-crash population size.

    Return the maximum and minimum pre-crash population sizes,
    and the corresponding time steps.
    """
    pre_crash_pop = popn.analytics_base.tumoursize[:treatmt.select_time]
    min_val = min(pre_crash_pop)
    max_val = max(pre_crash_pop)
    min_time = pre_crash_pop.index(min_val)
    max_time = pre_crash_pop.index(max_val)
    return min_val, min_time, max_val, max_time
