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



