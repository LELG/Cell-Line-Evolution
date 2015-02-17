#!/bin/bash

# if we're on the cluster, make sure correct PYTHONPATH is exported
if [[ ! -z $PBS_SERVER ]] && [[ $PBS_SERVER == "bioinf-head.petermac.org.au" ]]; then
  export PYTHONPATH=/usr/local/cluster/all_arch/python_libraries/production/lib/python2.7/site-packages
  export TERM=xterm
fi

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="config_file"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=1
if [ $# -ne $num_expected_args ]; then
  usage-exit "Incorrect number of args provided."
fi

# read parameters from config file
echo "================================================================================"
echo "Reading parameters from config file ..."
source $1
echo "Done"

# get appropriate padding for run_dir
# (this just finds the length of runs_per_param_set as a string)
runpadding=${#runs_per_param_set}

run_number=1
while [ $run_number -le $runs_per_param_set ]; do
  # note need to 'pad' run directories for proper ordering (001, 002 ... etc)
  run_dir=$(printf "%s/%0*d" $param_set_dir $runpadding $run_number)

  if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
  fi

  echo "------------------------------------------"
  echo
  echo "Starting simulation ..."
  echo
  echo "SIMULATION TEST GROUP: "$test_group
  echo "PARAMETER SET: "$param_set" of "$num_param_sets
  echo "RUN: "$run_number" of "$runs_per_param_set
  echo
  echo "------------------------------------------"
  echo

  # build up full list of parameters to main.py
  sim_params="--test_group $test_group --param_set $param_set --run_number $run_number"
  sim_params="$sim_params --test_group_dir $test_group_dir"
  sim_params="$sim_params --param_set_dir $param_set_dir"
  sim_params="$sim_params --run_dir $run_dir"
  sim_params="$sim_params --max_cycles $max_cycles --max_size_lim $max_size_lim"
  sim_params="$sim_params -p $proliferation_rate -d $death_rate -m $mutation_rate"
  sim_params="$sim_params --treatment_type $treatment_type"
  sim_params="$sim_params --decay_type $decay_type --decay_rate $decay_rate"
  sim_params="$sim_params --select_time $select_time --select_pressure $selective_pressure"
  sim_params="$sim_params --prob_mut_pos $prob_mut_pos --prob_mut_neg $prob_mut_neg"
  sim_params="$sim_params --prob_inc_mut $prob_inc_mut --prob_dec_mut $prob_dec_mut"
  sim_params="$sim_params --init_size $init_size"
  sim_params="$sim_params --init_diversity $init_diversity $sub_file"
  sim_params="$sim_params --scale $scale --mscale $mscale"
  sim_params="$sim_params $r_flag $m_flag $z_flag $np_flag"

  # run simulation with the full set of parameters
  python2.7 main.py $sim_params

  run_number=$((run_number+1))
done