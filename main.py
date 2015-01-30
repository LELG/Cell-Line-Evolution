""" 
LE MODIFICATIONS TO THE GA EXECUTION

Genomic Instability Simulation 

Simulate growth of and mutation of cancer cells 
and introduction of selective pressure

 
    python main.py -d 0.03 -p 0.04 -m 0.001 --loops 1000000 --maxsize_lim 1000000 --prolif_lim 0.0 --init_size 25 -f simple/1-1/ -s 0.01 -t 400000 --prob_mut_pos 0.01 --prob_mut_neg 0.99 --prob_inc_mut 0.0 --prob_dec_mut 0.0 --scale 0.5 --mscale 1.0 --M -n simple -g simple/1 --init_diversity 1 --sub_file /Users/Luis/Documents/MSc/Thesis/Code/src/AB/popln/heterogeneous_dat/zero_inc.dat --Z --NP
    testgroup='simple/1', testname='simple')
 

"""

from __future__ import print_function
import argparse
import population
import os


'''
    This function cast a dict into a Namespace
'''
class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)

def run_simulation(opt):
    """ 
    Create new population and analytics set to store results 
    Print info about population then begin simulation cycle 
    
    """
    
    population_base = population.Population(opt)
    #population_base.info()
    val = population_base.cycle(opt)
    print('LE ->Main ', val)
   
    #while not population_base.cycle(opt):
     #   population_base = population.Population(opt)
       # population_base.info()
        #print("restarting simulation - did not grow")

def ComputeFitness(FamilyMemeber):
    print ('Running simulation from '+ FamilyMemeber.g )
    return ( callFromGA(FamilyMemeber.genome) )

def callFromGA(Parameters_dict):
    path = os.getcwd() + '/heterogeneous_dat/zero_inc.dat'
    L = {'A':False, 
         'M':True, 
         'NP':True, 
         'R':False, 
         'Z':True, 
         'die':Parameters_dict['DeathRate'],
         'filename':'simple/1-1/',
         'init_diversity':Parameters_dict['InitialDiversity'],
         'init_size':Parameters_dict['InitSize'],
         'loops':Parameters_dict['loops'],
         'maxsize_lim':Parameters_dict['MaxSize'],
         'mscale':Parameters_dict['MutScale'],
         'mut':Parameters_dict['MutationRate'],
         'mut_amount_change':None,
         'mutagenic_pressure':0.0,
         'pro':Parameters_dict['ProliferationRate'],
         'prob_dec_mut':Parameters_dict['P_DecMutRate'],
         'prob_inc_mut':Parameters_dict['P_IncMutRate'],
         'prob_mut_neg':Parameters_dict['P_DelMutation'],
         'prob_mut_pos':Parameters_dict['P_BenMutation'],
         'prolif_lim':0.0,
         'scale':Parameters_dict['scale'],
         'select_pressure':Parameters_dict['SelectivePressure'],
         'select_time':400000,
         'sub_file':path,
         'testgroup':'simple/1',
         'testname':'simple'
         }
    opt = Bunch(L)
    print('Running with Parameters', Parameters_dict)
    FitnessValues = run(Parameters_dict, opt)
    return (FitnessValues)

    #opt = Bunch(Parameters_dict)
    
def run(Parameters_dict,opt):
    Fitnessvalues =  list()
    population_base = population.Population(opt)
    for i in xrange(Parameters_dict['Permutation']):
        val = population_base.cycle(opt) 
        if val == 0:
            val = 0.00000000001
        Fitnessvalues.append(val)

    print('Simulation values are ', Fitnessvalues)
    return (Fitnessvalues)



if __name__ == '__main__':
    callFromGA(dict())
"""
def main():
     Read simulation parameters from command line
    
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

    print(opt)

    #run_simulation(vars(opt))
    run_simulation(opt)

main()
"""
