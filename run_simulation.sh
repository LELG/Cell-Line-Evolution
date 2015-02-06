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

# if we're on the cluster, make sure correct PYTHONPATH is exported
if [[ ! -z $PBS_SERVER ]] && [[ $PBS_SERVER == "bioinf-head.petermac.org.au" ]]; then
  export PYTHONPATH=/usr/local/cluster/all_arch/python_libraries/production/lib/python2.7/site-packages
fi

# first, check for correct invocation
E_WRONG_ARGS=85
script_parameters="testgroup_name [config_file]"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=2
if [ $# -gt $num_expected_args ]; then
  usage-exit "Too many parameters provided"
elif [ $# -eq 0 ]; then
  usage-exit "Need test name as first argument"
fi

# check that supplied testname will make a valid directory name
# i.e., that it does not contain slash or backslash
E_INVALID_FNAME=60
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
  echo "Error: test group name cannot contain slash/backslash"
  exit $E_INVALID_FNAME
else
  testname=$1
fi

# make directory for today's date, unless it already exists
todaydir="results/$(date +'%Y-%m-%d')"

if [ ! -d "$todaydir" ]
then
  echo "Creating test directory for "$(date +'%Y-%m-%d')" ..."
  mkdir -p "$todaydir"
  echo "Done"
fi

# make main directory for this test group
maintestdir="$todaydir/$testname"
if [ ! -d $maintestdir ]
then
  echo "Creating main directory for test group "$testname" ..."
  mkdir -p $maintestdir
  echo "Done"
else
  echo "Warning: results for test group "$testname" already exist"
  i=1
  while [ -d $maintestdir-$i ] ; do
    let i++
  done
  maintestdir="$maintestdir-$i"
  echo "Creating directory test group "$maintestdir" ..."
  mkdir -p $maintestdir
  echo "Done"
fi

# read parameters from config file, if one has been specified
if [ -z $2 ]; then
  echo "No config file supplied, reading default parameters ..."
  source "default.conf"
else
  echo "Reading parameters from config file "$2" ..."
  source $2
  echo "Done"
fi

# make relevant directories
# run one / multiple sims

#pre_header="mutation_values: "$mutation_values", tests per group: "$runs_per_param_set
#echo $pre_header >> $testname/"results.dat"

# TODO check whether these files are still necessary
touch $maintestdir/"middropdata.dat"
touch $maintestdir/"enddropdata.dat"
#echo "TESTSIZE, etc"

# the following code assumes multiple mutation values and multiple tests
# TODO make this more general

param_set=1
for mutation_rate in $mutation_values; do
  run_number=1
  while [ $run_number -le $runs_per_param_set ]; do
    run_dir=$maintestdir/$param_set-$run_number/
    testgroup=$testname/$param_set

    if [ ! -d $run_dir ]; then
      mkdir -p $run_dir
    fi

    echo "Parameter set "$param_set"; run "$run_number" of "$runs_per_param_set

    python2.7 main.py --param_set $param_set --run_number $run_number -d $death_rate -p $proliferation_rate -m $mutation_rate --max_cycles $max_cycles --maxsize_lim $maxsize_lim --prolif_lim 0.0 --run_dir $run_dir --maintestdir $maintestdir -s $selective_pressure -t $select_time --prob_mut_pos $prob_mut_pos --prob_mut_neg $prob_mut_neg --prob_inc_mut $prob_inc_mut --prob_dec_mut $prob_dec_mut --scale $scale --mscale $mscale -n $testname -g $testgroup $diversity $r_flag $m_flag $z_flag $np_flag

    run_number=$((run_number+1))
  done

  param_set=$(($param_set+1))

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
