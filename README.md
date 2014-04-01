# popln\_mut
==========

A program to simulate a heterogeneous tumour population growth, mutation and strong selective pressure.

#### Summary

The simulation is run by a python program. Several bash and python scripts are then used to run multiple simulations and stitch the data together to plot the result of multiple runs.

#### Requirements

- Python 2.x with matplotlib, numpy
- Unix environment for shell scripts 

#### Outputs

The main program outputs
 - Two .xml files with an entry for each clone alive 
    - precrash
    - postcrash
 - Plots for single run for
    - precrash
    - postcrash
 - CSV Data of simulation

## Instructions

### 1. Gather all the data and scripts required.

 1. First Clone this repository
    
        git clone 

 2. Install matplotlib, numpy

### 2. Customise the `tests\parameter.sh` file

This file has all the options required for the simulation. One can duplicate the parameter file,
such as `parameter_mutation.sh` and change the values for a difference test. 
The test directory might end up filled with multiple tests.

The parameter file runs the `main.py` file with command line arguments specified within.
It uses one or more loops to iterate through a range of test values, given as command line arguments.

The simulation begins either with a single homogeneous clone or with a series of clones (heterogeneous
population). The heterogeneous population is read in from a hetpop.dat file, specified in parameters.sh.

The following parameters can be changed in the shell script:

        --filename 
        --testname
        --testgroup
        --loops - Number of Loops/Cycles in simulation
        --die - Death rate
        --pro - Proliferation rate
        --mut - Mutation rate
        --maxsize_lim - Maximum Tumour Size
        --prolif_lim  - Proliferation Limit #not used
        --init_size   - Initial Tumour Size
        --select_time - Selective Pressure Time
        --select_pressure - Amount of Selective Pressure
        --mutagenic_pressure - Mut drop/inc during selection event
        --scale        - scaling parameter, partly hard coded at the moment
        --prob_mut_pos - Probability that mutation with be beneficial
        --prob_mut_neg - Probability that mutation will be deleterious
        --mut_amount_change - Amount of max mutation change
        --prob_inc_mut - Probability of increasing mutation rate
        --prob_dec_mut - Probability of decreasing mutation rate
        --init_diversity - Initial Diversity of Population
        --R - Output data in R format instead of matplotlib
        --M - Introduce selective pressure automatically at max size
        --A - Broken
        --Z - Prune tree while running for larger executions
    

##### Homogeneous Population

See parametersHom.sh

##### Heterogeneous Population

See parametersHet.sh

### 3. Run the simulation(s)

Run the parameter file, specifying a unique test name

        ./parameter.sh TEST_NAME

The output of the file with n tests for each combination of parameters, x.

        ./TEST_NAME/1-1
        ./TEST_NAME/1-2
        ./TEST_NAME/1-3
        ./TEST_NAME/2-1
        ...
        ./TEST_NAME/x-n

There will also be a 

        ./TEST_NAME/1-1/_phylo.xml

Note, I will refer to a run as a single instance of `main.py` while a test involves running `parameter.sh` with potentially many runs.

The plots produced by a single run are:

* Populations - Tumour size against time (number of discrete cycles)
* Subpopulations - Plot of the number of clones (alive) against time. 
* Mutation Rate - Plot of average mutation rate against time
* Mutation Distribution – Mid and End - Histogram of the mutation rates distributed acoss the total number of cells. Two graphs are plotted together, one showing the distribution before the crash, and the second showing the distribution of mutation rates after the crash. Thus we can see the change and trends in mutation rate.
* Proliferation Rate - Plot of the effective proliferation rate, that is the average of the proliferation rate for each cell, taking into account a reduction in proliferation due to population constraints as well as due to selective pressure.
* Proliferation Distribution – Mid and End - Histogram of the proliferation rates distributed acoss the total number of cells. Two graphs are plotted together, one showing the distribution before the crash, and the second showing the distribution of mutation rates after the crash. Thus we can see the change and trends in proliferation rate.
* Cell Line Graph - Each vertical line indicates the emergence of a clone. When the line ends, that clone has died out (or maybe it to the end of the simulation). Graph is sorted by time of initial emergence of that clone.
* Alleles - Shows the distribution of alleles (or mutations) across the population. That is, if one mutation is found in a lot of desendents, it will show up as one of the larger allele frequencies.



### 4. Compile results from multiple runs

There are several options for compiling results. The script `compact.py` reads the combined output from the results.dat file found in the top level of a test. The results are for all the runs under one `parameter.sh` file.

    ./TEST_NAME/results.dat

Running `compact.py` groups this, so instead of x * n rows, there are only x rows, which summarise n runs each. It also has lateoutput format `--latex` and a headerless `--plot` argument. Normal execution is:

    compact.py --filename FILENAME 

The other way of summarising data across heterogenous populations for one test.

### 5. Plot results from multiple runs - Homogenous

    ./summary_winner.sh TEST_NAME

### 5. Plot results from multiple runs - Heterogeneous

Create a summary of heterogeneous run, currently hard coded to retrieve the results of each individual colour. Use summary script to retrieve values.
    ./summary_xml_zero.sh TEST_NAME

Make a bar graph.

    ./pr_bar_manual.py



