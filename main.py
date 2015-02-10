"""
Genomic Instability Simulation

Simulate growth of and mutation of cancer cells
and introduction of selective pressure

"""
from __future__ import print_function

# Force matplotlib not to use any Xwindows backend
# to avoid error on the PMCI cluster. Note that
# this import and call to matplotlib.use() must
# appear before any other matplotlib import
import matplotlib
matplotlib.use('Agg')

import argparse
import os
import csv
import simulator
# import population


def main():
    """
    Parse simulation parameters, create results file, and run simulation.

    This function gets simulation parameters from the
    command line; creates two new summary files (one for storing
    all results from this group of tests, the other containing
    only results from a single parameter set); and starts
    the simulation with the given parameter set. The parameters
    are as follows:

    Returns
    -------
    None.

    See also
    --------
    parse_cmd_line_args : describes command line arguments
    """
    opt = parse_cmd_line_args()

    # calculate proliferation limit and store in parameter set
    opt.prolif_lim = opt.pro - opt.die

    # create results files, if they don't already exist
    tgroup_summ_path = "{0}/{1}_results.csv".format(opt.test_group_dir,
                                                    opt.test_group)
    pset_summ_path = "{0}/{1}/{2}_{1}_results.csv".format(opt.test_group_dir,
                                                          opt.param_set,
                                                          opt.test_group)
    if not os.path.exists(tgroup_summ_path):
        initialise_results(tgroup_summ_path)
    if not os.path.exists(pset_summ_path):
        initialise_results(pset_summ_path)

    # finally, run simulation
    sim = simulator.Simulator(opt)
    sim.popn.info()
    sim.run()


def parse_cmd_line_args():
    """
    Parse simulation parameters.

    This function reads in simulation parameters from
    command line arguments, storing them in a namespace object.
    The parameters are as follows:

    test_group : string
        Identifier for current group of simulations.
    test_group_dir : string
        Main directory for this group of simulations.
    param_set : string
        Identifier of a given parameter set.
    run_number : integer
        Counter that tracks how many times this
        particular parameter set has been run in
        this set of simulations.

    max_cycles : integer
        Maximum number of time steps in simulation
    max_size_lim : integer
        Maximum tumour size. This is the tumour size
        which triggers, first, selective pressure, and
        second, the end of the simulation.

    pro : float
        Initial proliferation rate (homogeneous population)
    die : float
        Initial death rate
    mut : float
        Initial mutation rate (homogeneous population)
    init_size : integer
        Initial tumour size (homogeneous population)

    select_time : integer
        Time step to introduce selective pressure,
        if it has not already been automatically triggered.
    select_pressure : float
        Amount of Selective Pressure
    mutagenic_pressure : float
        Change in mutation rate during selection event

    prob_mut_pos : float
        Probability that mutation with be beneficial
    prob_mut_neg : float
        Probability that mutation will be deleterious
    prob_inc_mut : float
        Probability of increasing mutation rate
    prob_dec_mut : float
        Probability of decreasing mutation rate

    init_diversity : bool
        Initial diversity of population
        (homogeneous or heterogeneous)
    sub_file : string
        Name of file containing list of
        subpopulations for an initially
        heterogeneous population.

    scale : float
        Scaling parameter, partly hard coded at the moment
    mscale : float
        Mutation rate scaling parameter, partly hard coded at the moment

    --R : bool
        Output data in R format instead of matplotlib
    --M : bool
        Introduce selective pressure automatically at max size
    --Z : bool
        Prune tree while running for larger executions
    --NP: bool
        Don't print plots for this simulation

    Returns
    -------
    opt : an argparse.Namespace() object containing
        the values of the above parameters.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--test_group', default='default_test_group')
    parser.add_argument('--test_group_dir')
    parser.add_argument('--param_set', type=str, default='1')
    parser.add_argument('--run_number', type=int, default=1)

    parser.add_argument('-l', '--max_cycles', type=int)
    parser.add_argument('-x', '--max_size_lim', type=int)

    parser.add_argument('-p', '--pro', type=float)
    parser.add_argument('-d', '--die', type=float)
    parser.add_argument('-m', '--mut', type=float)

    parser.add_argument('-t', '--select_time', type=int, default=400000)
    parser.add_argument('-s', '--select_pressure', type=float, default=0.01)
    parser.add_argument('-u', '--mutagenic_pressure', type=float, default=0.0)

    parser.add_argument('--prob_mut_pos', type=float)
    parser.add_argument('--prob_mut_neg', type=float)
    parser.add_argument('--prob_inc_mut', type=float)
    parser.add_argument('--prob_dec_mut', type=float)

    parser.add_argument('-z', '--init_size', type=int, default=25)
    parser.add_argument('--init_diversity', type=int, default=0)
    parser.add_argument('--sub_file', default='zero_inc.csv')

    parser.add_argument('-c', '--scale', type=float, default=1000.0)
    parser.add_argument('-e', '--mscale', type=float, default=1000.0)

    parser.add_argument('--R', action="store_true", default=False)
    parser.add_argument('--M', action="store_true", default=False)
    parser.add_argument('--Z', action="store_true", default=False)
    parser.add_argument('--NP', action="store_true", default=False)

    opt = parser.parse_args()
    return opt


def initialise_results(filepath):
    """
    Initialise results file for this simulation.

    Attempt to create a new summary file in the
    specified directory, raising an exception if
    the directory does not exist.

    This function specifies the column structure of the
    results file; at some later date, it might make
    sense to specify the structure elsewhere, and pass
    that as an argument to this function.

    Parameters
    ----------
    filepath : string
        Location of the results file to be created.

    Returns
    -------
    None. Writes to file.

    """

    columns = ('param_set', 'run_number',
               'went_through_crash',
               'recovered', 'recov_type', 'recov_percent',
               'prolif_rate', 'death_rate', 'mut_rate',
               'select_time', 'select_pressure',
               'prob_ben_mut', 'prob_del_mut',
               'prob_mut_incr', 'prob_mut_decr',
               'pop_size', 'num_clones',
               'avg_mut_rate_at_end', 'avg_prolif_rate_at_end',
               'elapsed_time', 'elapsed_cycles',
               'pre_crash_min', 'pre_crash_min_time',
               'pre_crash_max', 'pre_crash_max_time',
               'post_crash_min', 'post_crash_min_time',
               'post_crash_max', 'post_crash_max_time')

    try:
        results_file = open(filepath, "w")
    except IOError:
        raise

    results_writer = csv.writer(results_file)
    results_writer.writerow(columns)
    results_file.close()

'''
def run_simulation(opt):
    """
    Run simulation.

    Create a new population (including associated
    analytics object); print details of population
    parameters to stdout; then begin simulation cycle.

    If the population dies out prematurely, restart
    the simulation.

    Parameters
    ----------
    opt : argparse.Namespace object
        Set of simulation parameters, originally
        specified as command line arguments

    Returns
    -------
    No return value. Creates several output files.

    See Also
    --------
    population and subpopulation, for details
    on output files.
    """

    population_base = population.Population(opt)
    population_base.info()

    while not population_base.cycle(opt):
        population_base = population.Population(opt)
        population_base.info()
        print("restarting simulation - did not grow")
'''

if __name__ == '__main__':
    main()
