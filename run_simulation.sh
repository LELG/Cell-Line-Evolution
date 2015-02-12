#!/bin/bash
#
# MAIN SCRIPT FOR RUNNING TUMOUR SIMULATION
#

# if we're on the cluster, make sure correct PYTHONPATH is exported
if [[ ! -z $PBS_SERVER ]] && [[ $PBS_SERVER == "bioinf-head.petermac.org.au" ]]; then
  export PYTHONPATH=/usr/local/cluster/all_arch/python_libraries/production/lib/python2.7/site-packages
  export TERM=xterm
fi

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="test_group_name [config_file]"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=2
if [ $# -gt $num_expected_args ]; then
  usage-exit "Too many parameters provided"
elif [ $# -eq 0 ]; then
  usage-exit "Need test group name as first argument"
fi

# check that supplied testname will make a valid directory name
# i.e., that it does not contain slash or backslash
E_INVALID_FNAME=60
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
  echo "Error: test group name cannot contain slash/backslash"
  exit $E_INVALID_FNAME
else
  test_group=$1
fi

# make directory for today's date, unless it already exists
today=$(date +'%Y-%m-%d')
today_dir="results/$today"

if [ ! -d "$today_dir" ]
then
  echo "Creating results directory for "$today" ..."
  mkdir -p "$today_dir"
  echo "Done"
fi

# make main directory for this test group
test_group_dir="$today_dir/$test_group"
if [ ! -d $test_group_dir ]
then
  echo "Creating main directory for test group "$test_group" ..."
  mkdir -p $test_group_dir
  echo "Done"
else
  echo "Warning: results for test group "$test_group" already exist"
  i=1
  while [ -d "$test_group_dir($i)" ] ; do
    let i++
  done
  test_group_dir="$test_group_dir($i)"
  echo "Creating test group directory: "$test_group_dir" ..."
  mkdir -p $test_group_dir
  echo "Done"
fi

# read parameters from config file, if one has been specified
# TODO sort out which params come from .conf, which from PBS script
if [ -z $2 ]; then
  echo "No config file supplied, reading default parameters ..."
  source "default.conf"
else
  echo "Reading parameters from config file "$2" ..."
  source $2
  echo "Done"
  if [ ! -f $test_group_dir/$2 ]; then
    echo "Copying config file to results directory ..."
    cp $2 $test_group_dir
    echo "Done"
  fi
fi

# get appropriate padding for run_dir
# (this just finds the length of runs_per_param_set as a string)
runpadding=${#runs_per_param_set}

# TODO check whether these files are still necessary
touch $test_group_dir/"middropdata.csv"
touch $test_group_dir/"enddropdata.csv"

# the following code assumes multiple mutation values and multiple tests
# TODO delegate this looping to PBS script

param_set=1
for mutation_rate in $mutation_values; do
  param_set_dir="$test_group_dir/$param_set"

  if [ ! -d $param_set_dir ]; then
    mkdir -p $param_set_dir
  fi

  run_number=1
  while [ $run_number -le $runs_per_param_set ]; do
    # note need to 'pad' run directories for proper ordering (001, 002 ... etc)
    run_dir=$(printf "%s/%0*d" $param_set_dir $runpadding $run_number)

    if [ ! -d $run_dir ]; then
      mkdir -p $run_dir
    fi

    echo "TEST GROUP: "$test_group
    echo "Parameter set: "$param_set
    echo "Run "$run_number" of "$runs_per_param_set

    # build up full list of parameters to main.py
    sim_params="--test_group $test_group --param_set $param_set --run_number $run_number"
    sim_params="$sim_params --test_group_dir $test_group_dir"
    sim_params="$sim_params --param_set_dir $param_set_dir"
    sim_params="$sim_params --run_dir $run_dir"
    sim_params="$sim_params --max_cycles $max_cycles --max_size_lim $max_size_lim"
    sim_params="$sim_params -p $proliferation_rate -d $death_rate -m $mutation_rate"
    sim_params="$sim_params --select_time $select_time --select_pressure $selective_pressure"
    sim_params="$sim_params --prob_mut_pos $prob_mut_pos --prob_mut_neg $prob_mut_neg"
    sim_params="$sim_params --prob_inc_mut $prob_inc_mut --prob_dec_mut $prob_dec_mut"
    sim_params="$sim_params $init_size"
    sim_params="$sim_params $init_diversity $sub_file"
    sim_params="$sim_params --scale $scale --mscale $mscale"
    sim_params="$sim_params $r_flag $m_flag $z_flag $np_flag"

    # run simulation with the full set of parameters
    python2.7 main.py $sim_params

    run_number=$((run_number+1))
  done

  param_set=$(($param_set+1))

  #python plot_circles.py -f $param_set_dir"_circles.dat"
  #python plot_circles.py -f $param_set_dir"_circles_all.dat"
done

#python compact.py -f $testname/results.dat > $testname/summary.dat
#echo "about to run test"
#echo python -c "import dropdata; dropdata.read_drop("$testname")"
#python -c "import dropdata; dropdata.read_drop("\"$testname"\")"
#echo $testgroup
