"""
Parse command line arguments, and start the simulation.

Authors
-------
Andrew Bakshi : andrew.bakshi@gmail.com
Yoshua Wakeham : yoshwakeham@gmail.com

Date created
------------
01 April 2014
"""
from __future__ import print_function

# Force matplotlib not to use any Xwindows backend,
# to avoid an error on the PMCI cluster. Note that
# this import and call to matplotlib.use() must
# appear before any other matplotlib import
import matplotlib
matplotlib.use('Agg')

import argparse
import os
import csv
import simulator


def main():
    """
    Create results files and run simulation with given params.

    Get simulation parameters from the command line;
    create summary files (one for storing all results
    from this group of tests, the other only results
    from this parameter set); run the simulation with
    the given parameters.

    Args
    ----
    None

    Returns
    -------
    None
    """
    # get parameters
    opt = parse_cmd_line_args()

    # TODO move this initialisation to simulator?
    initialise_results(opt)

    # create and run simulation
    sim = simulator.Simulator(opt)
    sim.print_info()
    sim.run()


def parse_cmd_line_args():
    """
    Parse simulation parameters from command line.

    Read in simulation parameters from
    command line arguments, storing them in a
    namespace object.

    Args
    ----
    None. The accepted command line arguments
    are as follows:

    test_group : string
        Identifier for current group of simulations.
    param_set : string
        Identifier of a given parameter set.
    run_number : integer
        Counter that tracks how many times this
        particular parameter set has been run in
        this set of simulations.

    test_group_dir : string
        Main directory for this group of simulations.
    param_set_dir : string
        Directory for this particular parameter set.
    run_dir : string
        Directory for this specific replicate run.

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

    treatment_type : string
        Treatment regime (e.g. single dose, metronomic, adaptive)
    decay_type : string
        Function for specifying shape of selective pressure
        decay (e.g. constant, linear, exponential)
    decay_rate : float
        Constant controlling speed of treatment decay
    treatment_freq : int
        Number of cycles between doses for (regular) multi-dose treatments
    adaptive_increment : float
        Amount to increase/decrease dosage by in adaptive treatment
    adaptive_threshold : float
        Percentage change in tumour size that will trigger change
        in adaptive treatment dose
    select_time : integer
        Time step to introduce treatment, if it has
        not already been automatically triggered
    select_pressure : float
        Initial quantity of selective pressure
        introduced by treatment
    mutagenic_pressure : float
        Change in mutation rate during selection event

    NOTE: As currently implemented, mutagenic_pressure is
        considered to be a multiplicative factor, i.e. clones
        with a higher mutation rate experience higher mutagenic
        pressure. In future this may change to be an additive
        pressure, analogous to proliferation pressure.

    --resistance : bool
        Allow for pre-existing resistance mutations
    num_resist_mutns : int
        Deterministically specify number of resistance mutations.
        If not specified, number of mutations will be determined
        stochastically as function of tumour size and number of
        existing deleterious/neutral mutations.

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
    An argparse.Namespace() containing the
    supplied values of the above parameters.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--test_group')
    parser.add_argument('--param_set')
    parser.add_argument('--run_number', type=int)

    parser.add_argument('--test_group_dir')
    parser.add_argument('--param_set_dir')
    parser.add_argument('--run_dir')

    parser.add_argument('-l', '--max_cycles', type=int)
    parser.add_argument('-x', '--max_size_lim', type=int)

    parser.add_argument('-p', '--pro', type=float)
    parser.add_argument('-d', '--die', type=float)
    parser.add_argument('-m', '--mut', type=float)

    parser.add_argument('--treatment_type', default='single_dose')
    parser.add_argument('--decay_type', default='constant')
    parser.add_argument('--decay_rate', type=float, default=0.0)
    parser.add_argument('--treatment_freq', type=int, default=100)
    parser.add_argument('--adaptive_increment', type=float, default=0.001)
    parser.add_argument('--adaptive_threshold', type=float, default=0.025)
    parser.add_argument('-t', '--select_time', type=int)
    parser.add_argument('-s', '--select_pressure', type=float)
    parser.add_argument('-u', '--mutagenic_pressure', type=float, default=0.0)

    parser.add_argument('--resistance', action="store_true", default=False)
    parser.add_argument('--num_resist_mutns', type=int, default=-1)

    parser.add_argument('--prob_mut_pos', type=float)
    parser.add_argument('--prob_mut_neg', type=float)
    parser.add_argument('--prob_inc_mut', type=float)
    parser.add_argument('--prob_dec_mut', type=float)

    parser.add_argument('-z', '--init_size', type=int, default=25)
    parser.add_argument('--init_diversity', type=int)
    parser.add_argument('--sub_file')

    parser.add_argument('-c', '--scale', type=float)
    parser.add_argument('-e', '--mscale', type=float)

    parser.add_argument('--r_output', '--R',
                        action="store_true", default=False)
    parser.add_argument('--auto_treatment', '--M',
                        action="store_true", default=False)
    parser.add_argument('--prune_clones', '--Z',
                        action="store_true", default=False)
    parser.add_argument('--no_plots', '--NP',
                        action="store_true", default=False)

    return parser.parse_args()


def initialise_results(opt):
    """Create summary files, unless they already exist."""
    tgroup_summary_path = "{0}/{1}_results.csv".format(opt.test_group_dir,
                                                       opt.test_group)
    pset_summary_path = "{0}/{1}_{2}_results.csv".format(opt.param_set_dir,
                                                         opt.test_group,
                                                         opt.param_set)
    if not os.path.exists(tgroup_summary_path):
        create_results_file(tgroup_summary_path)
    if not os.path.exists(pset_summary_path):
        create_results_file(pset_summary_path)


def create_results_file(filepath):
    """
    Create a results file.

    Attempt to create a new simulation summary
    file in the specified directory, raising an
    exception if the directory does not exist.

    Note
    ----
    This function specifies the column structure of the
    results file; at some later date, it might make
    sense to specify the structure elsewhere, and pass
    that as an argument to this function.

    Args
    ----
    filepath : string
        Location of the results file to be created.

    Returns
    -------
    None

    Raises
    ------
    IOError: If `filepath` is invalid.
    """

    columns = ('param_set', 'run_number',
               'went_through_crash',
               'recovered', 'recov_type', 'recov_percent',
               'prolif_rate', 'death_rate', 'mut_rate',
               'crash_time', 'select_pressure',
               'max_select_pressure',
               'prob_ben_mut', 'prob_del_mut',
               'prob_mut_incr', 'prob_mut_decr',
               'pop_size', 'num_clones',
               'dom_clone_proportion', 'avg_depth',
               'avg_mut_rate_at_end', 'avg_prolif_rate_at_end',
               'elapsed_time', 'elapsed_cycles',
               'total_mutations',
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


if __name__ == '__main__':
    main()
