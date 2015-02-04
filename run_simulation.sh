#!/bin/bash

# MAIN SCRIPT FOR RUNNING TUMOUR SIMULATION
#
# This script gets the test parameters from
# a config file, and supplies them to the main
# Python routine. It co-ordinates the running of
# multiple simulations; and hopefully in future
# will co-ordinate the complicated qsub scheduling
# that may be required for running 100s or 1000s of
# tests.


# first, check for correct invocation
E_WRONG_ARGS=85
num_expected_args=2
script_parameters="testgroup_name [config_file]"

if [ $# -gt $num_expected_args ]; then
  echo "Usage: ./`basename $0` $script_parameters"
  echo "Too many parameters provided"
  exit $E_WRONG_ARGS
elif [ $# -eq 0 ]; then
  echo "Usage: ./`basename $0` $script_parameters"
  echo "Need test name as first argument"
  exit $E_WRONG_ARGS
fi

# check that supplied testname will make a valid directory name
E_INVALID_FNAME=60
if [ "${1/'/'}" != "$1" ] ; then
  echo "Invalid test group name: $1"
  exit $E_INVALID_FNAME
else
  testname=$1
fi

# make main directory for this test group
if [ ! -d $testname ]
then
  echo "Creating main test directory ..."
  mkdir $testname
  echo "Done"
fi

# read parameters from config file, if one has been specified
if [ -z $2 ]; then
  echo "No config file supplied, reading default parameters ..."
  source "default.conf"
else
  echo "Reading parameters from config file"$2"..."
  source $2
fi

# make relevant directories
# run one / multiple sims

#pre_header="mutation_values: "$mutation_values", tests per group: "$runs_per_param_set
#echo $pre_header >> $testname/"results.dat"

# TODO check whether these files are still necessary
touch $testname/"middropdata.dat"
touch $testname/"enddropdata.dat"
#echo "TESTSIZE, etc"

# the following code assumes multiple mutation values and multiple tests
# TODO make this more general

param_set=1
for mutation_rate in $mutation_values; do
  run_number=1
  while [ $run_number -le $runs_per_param_set ]; do
    filepath=$testname/$param_set-$run_number/
    testgroup=$testname/$param_set

    if [ ! -d $filepath ]; then
      mkdir $filepath
    fi

    echo "Parameter set "$param_set", run "$run_number" of "$runs_per_param_set

    python main.py --param_set $param_set --run_number $run_number -d $death_rate -p $proliferation_rate -m $mutation_rate --loops $loops --maxsize_lim $maxsize_lim --prolif_lim 0.0 --init_size $initial_size -f $filepath -s $selective_pressure -t $select_time --prob_mut_pos $prob_mut_pos --prob_mut_neg $prob_mut_neg --prob_inc_mut $prob_inc_mut --prob_dec_mut $prob_dec_mut --scale $scale --mscale $mscale -n $testname -g $testgroup --init_diversity $initial_diversity $sub_file $r_flag $m_flag $z_flag $np_flag

    run_number=$((run_number+1))
  done

  count=$(($count+1))

  #INDICATE AN IDENTICAL SET OF TESTS FINISHED
  #echo "END_GROUP" >> $testname/"results.dat"
  #python plot_circles.py -f $testgroup"_circles.dat"
  #python plot_circles.py -f $testgroup"_circles_all.dat"
done

#python compact.py -f $testname/results.dat > $testname/summary.dat
#echo "about to run test"
#echo python -c "import dropdata; dropdata.read_drop("$testname")"
#python -c "import dropdata; dropdata.read_drop("\"$testname"\")"
#echo $testgroup
