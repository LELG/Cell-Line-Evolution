""" Genomic Instability Simulation

Simulate growth of and mutation of cancer cells 
and introduction of selective pressure

"""

from __future__ import print_function
import argparse
import os
import errno
import population

def run_simulation(opt):
    """ 
    Create new population and analytics set to store results 
    Print info about population then begin simulation cycle 
    
    """
    
    population_base = population.Population(opt)
    population_base.info()
   
    while not population_base.cycle(opt):
        population_base = population.Population(opt)
        population_base.info()
        print("restarting simulation - did not grow")


def make_path_unless_exists(path):
    """Make a directory, unless it already exists"""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def main():
    """ Read simulation parameters from command line
    
    filename 
    testname
    testgroup
    loops - Number of Loops/Cycles in simulation
    die - Death rate
    pro - Proliferation rate
    mut - Mutation rate
    maxsize_lim - Maximum Tumour Size
    prolif_lim  - Proliferation Limit #not used
    init_size   - Initial Tumour Size
    select_time - Selective Pressure Time
    select_pressure - Amount of Selective Pressure
    mutagenic_pressure - Mut drop/inc during selection event
    scale        - scaling parameter, partly hard coded at the moment
    prob_mut_pos - Probability that mutation with be beneficial
    prob_mut_neg - Probability that mutation will be deleterious
    mut_amount_change - Amount of max mutation change
    prob_inc_mut - Probability of increasing mutation rate
    prob_dec_mut - Probability of decreasing mutation rate
    init_diversity - Initial Diversity of Population
    --R - Output data in R format instead of matplotlib
    --M - Introduce selective pressure automatically at max size
    --A - Broken
    --Z - Prune tree while running for larger executions
    
    """

    params = ['>', 'filename', 'went through crash', 'recovered',
              'recovery type', 'recovery percent',
              'pro', 'die', 'mut', 'select time',
              'select pressyre', 'prob pos mut', 'prob neg mut',
              'prob of increase in mut rate', 'prob of mut incr',
              'prob of mut decr', 'pop size', 'no. of clones',
              'avg mut rate at end', 'avg pro rate at end', 'time',
              'tumoursize', 'cycles', 'pre crash in', 'pre crash min time',
              'pre crash max', 'ax time', 'post crash min', 'min time',
              'post crash max', 'max time',
              'scale', 'mutation scale', 'size at crash']

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default='filename')
    parser.add_argument('--sub_file', default='sub_file.dat')
    parser.add_argument('-n', '--testname', default='testname')
    parser.add_argument('-g', '--testgroup', default='testname')
    parser.add_argument('-l', '--loops', type=int)
    parser.add_argument('-d', '--die', type=float)
    parser.add_argument('-p', '--pro', type=float)
    parser.add_argument('-m', '--mut', type=float)
    parser.add_argument('-x', '--maxsize_lim', type=int)
    parser.add_argument('-i', '--prolif_lim', type=float)
    parser.add_argument('-z', '--init_size', type=int)
    parser.add_argument('--init_diversity', type=int, default=0)
    parser.add_argument('-t', '--select_time', type=int)
    parser.add_argument('-s', '--select_pressure', type=float)
    parser.add_argument('-u', '--mutagenic_pressure', type=float, default=0.0)
    parser.add_argument('-c', '--scale', type=float, default=1000.0)
    parser.add_argument('-e', '--mscale', type=float, default=1000.0)
    parser.add_argument('--prob_mut_pos', type=float)
    parser.add_argument('--prob_mut_neg', type=float)
    parser.add_argument('--mut_amount_change', type=float)
    parser.add_argument('--prob_inc_mut', type=float)
    parser.add_argument('--prob_dec_mut', type=float)
    parser.add_argument('--R', action="store_true", default=False)
    parser.add_argument('--M', action="store_true", default=False)
    parser.add_argument('--A', action="store_true", default=False)
    parser.add_argument('--Z', action="store_true", default=False)
    parser.add_argument('--NP', action="store_true", default=False)
    opt = parser.parse_args()

    #run_simulation(vars(opt))
    run_simulation(opt)

main()

