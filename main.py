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
    sim.run(sim.start_cycle)


def parse_cmd_line_args():
    """
    Parse simulation parameters from command line.

    Read in simulation parameters from
    command line arguments, storing them in a
    namespace object.

    Args
    ----
    The accepted command line arguments
    are as follows:

    IDENTIFIERS
    ===========
    test_group : string
        Identifier for current group of simulations.
    param_set : string
        Identifier of a given parameter set.
    run_number : integer
        Counter that tracks how many times this
        particular parameter set has been run in
        this set of simulations.

    DIRECTORIES
    ===========
    test_group_dir : string
        Main directory for this group of simulations.
    param_set_dir : string
        Directory for this particular parameter set.
    run_dir : string
        Directory for this specific replicate run.

    BOUNDS
    ======
    max_cycles : integer
        Maximum number of time steps in simulation
    max_size_lim : integer
        Maximum tumour size. This is the tumour size
        which triggers, first, selective pressure, and
        second, the end of the simulation.

    TUMOUR PARAMETERS
    =================
    pro : float
        Initial proliferation rate (homogeneous population)
    die : float
        Initial death rate
    mut : float
        Initial mutation rate (homogeneous population)
    init_size : integer
        Initial tumour size (homogeneous population)
    init_diversity : bool
        Initial diversity of population
        (homogeneous or heterogeneous)
    sub_file : string
        Name of file containing list of
        subpopulations for an initially
        heterogeneous population.

    TREATMENT PARAMETERS
    ====================
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

    RESISTANCE PARAMETERS
    =====================
    --resistance : bool
        Allow for pre-existing resistance mutations
    num_resist_mutns : int
        Deterministically specify number of resistance mutations.
        If not specified, number of mutations will be determined
        stochastically as function of tumour size and number of
        existing deleterious/neutral mutations.
    resist_strength : float
        Determine the strength of the immunity to selective
        pressure conferred by a resistance mutation. This must be
        a float between 0 and 1.

    PROBABILITIES
    =============
    prob_mut_pos : float
        Probability that mutation with be beneficial
    prob_mut_neg : float
        Probability that mutation will be deleterious
    prob_inc_mut : float
        Probability of increasing mutation rate
    prob_dec_mut : float
        Probability of decreasing mutation rate

    SAVING/LOADING
    ==============
    * NOTE: the following options are mutually exclusive

    save_snapshot : bool
        Save a snapshot of the simulation
        population just before treatment
        is introduced.
    load_snapshot SNAPSHOT_ARCHIVE : string
        Load a population from a stored snapshot.
        SNAPSHOT_ARCHIVE specifies the filepath of the
        snapshot archive.

    SCALING
    =======
    scale : float
        Scaling parameter, partly hard coded at the moment
    mscale : float
        Mutation rate scaling parameter, partly hard coded at the moment

    MISCELLANEOUS
    =============
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

    # define a custom argparse action to enable cmd line options like
    #     --load_snapshot SNAPSHOT_ARCHIVE
    # rather than requiring the messier
    #     --load_snapshot --snapshot_archive SNAPSHOT_ARCHIVE
    # This code is based on the accepted answer at
    # http://stackoverflow.com/questions/8632354
    def store_flag_and_vars(varnames):
        """varnames is a list of variables to set based on cmd line input"""
        class customAction(argparse.Action):
            def __call__(self, parser, opt, values, option_string=None):
                setattr(opt, self.dest, True)
                if isinstance(values, list):
                    if len(values) != len(varnames):
                        errmsg = "expected {} varnames, got {}"
                        raise ValueError(errmsg.format(len(values), len(varnames)))
                    for i in range(len(values)):
                        setattr(opt, varnames[i], values[i])
                else:
                    if len(varnames) != 1:
                        errmsg = "expected one varname, got {}"
                        raise ValueError(errmsg.format(len(varnames)))
                    setattr(opt, varnames[0], values)
        return customAction

    parser = argparse.ArgumentParser()

    identifiers = parser.add_argument_group("identifiers")
    identifiers.add_argument('--test_group', required=True)
    identifiers.add_argument('--param_set', required=True)
    identifiers.add_argument('--run_number', required=True, type=int)

    directories = parser.add_argument_group("directories")
    directories.add_argument('--test_group_dir', required=True)
    directories.add_argument('--param_set_dir', required=True)
    directories.add_argument('--run_dir', required=True)

    bounds = parser.add_argument_group("bounds")
    bounds.add_argument('--max_cycles', type=int, default=int(1e5))
    bounds.add_argument('--max_size_lim', type=int, default=int(1e5))

    tumour_params = parser.add_argument_group("tumour params")
    tumour_params.add_argument('--pro', type=float, default=0.04)
    tumour_params.add_argument('--die', type=float, default=0.03)
    tumour_params.add_argument('--mut', type=float, default=0.001)
    mutex_tumour = tumour_params.add_mutually_exclusive_group(required=True)
    mutex_tumour.add_argument('--init_size', type=int, default=25)
    mutex_tumour.add_argument('--init_diversity',
                              action=store_flag_and_vars(['sub_file']),
                              default=False, metavar='SUB_FILE')

    treatmt_params = parser.add_argument_group("treatment params")
    treatmt_params.add_argument('--treatment_type', default='single_dose')
    treatmt_params.add_argument('--decay_type', default='constant')
    treatmt_params.add_argument('--decay_rate', type=float, default=0.0)
    treatmt_params.add_argument('--treatment_freq', type=int, default=100)
    treatmt_params.add_argument('--adaptive_increment', type=float, default=0.001)
    treatmt_params.add_argument('--adaptive_threshold', type=float, default=0.025)
    treatmt_params.add_argument('--select_time', type=int, default=400000)
    treatmt_params.add_argument('--select_pressure', type=float, default=0.01)
    treatmt_params.add_argument('--mutagenic_pressure', type=float, default=0.0)

    resist_params = parser.add_argument_group("resistance params")
    resist_params.add_argument('--resistance', action="store_true", default=False)
    resist_params.add_argument('--num_resist_mutns', type=int, default=-1)
    resist_params.add_argument('--resist_strength', type=float, default=1.0)

    probabilities = parser.add_argument_group("probabilities")
    probabilities.add_argument('--prob_mut_pos', type=float, default=0.01)
    probabilities.add_argument('--prob_mut_neg', type=float, default=0.99)
    probabilities.add_argument('--prob_inc_mut', type=float, default=0.0)
    probabilities.add_argument('--prob_dec_mut', type=float, default=0.0)

    saving_loading = parser.add_argument_group("saving/loading",
                                               "(mutually exclusive)")
    mutex_saveload = saving_loading.add_mutually_exclusive_group()
    mutex_saveload.add_argument('--save_snapshot',
                                action="store_true", default=False)
    mutex_saveload.add_argument('--load_snapshot',
                                action=store_flag_and_vars(['snapshot_archive']),
                                default=False, metavar="SNAPSHOT_ARCHIVE")

    scaling = parser.add_argument_group("scaling")
    scaling.add_argument('--scale', type=float, default=0.5)
    scaling.add_argument('--mscale', type=float, default=1.0)

    misc = parser.add_argument_group("miscellaneous")
    misc.add_argument('--r_output', '--R', action="store_true", default=False)
    misc.add_argument('--auto_treatment', '--M', action="store_true", default=False)
    misc.add_argument('--prune_clones', '--Z', action="store_true", default=False)
    misc.add_argument('--no_plots', '--NP', action="store_true", default=False)

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
               'generated_resist_mutns', 'surviving_resist_mutns',
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
