#!/bin/bash

# RUN GROUP OF SIMULATIONS
# where parameter sets are specified line-by-line in a space-separated config file

# check for correct invocation
E_WRONG_ARGS=85
script_parameters="test_group_name runs_per_param_set config_file"

function usage-exit () {
  echo "Usage: ./`basename $0` $script_parameters"
  echo >&2 $@
  exit $E_WRONG_ARGS
}

num_expected_args=3
if [ $# -ne $num_expected_args ]; then
  usage-exit "Incorrect number of parameters provided"
fi

if $(echo $2 | grep -E -q '^[0-9]+$'); then
  runs_per_param_set=$2
else
  usage-exit "Numeric arg required, got "$2" instead"
fi

if [ ! -f $3 ]; then
  usage-exit "Could not read config file: "$3
fi

# check that supplied testname will make a valid directory name
# i.e., that it does not contain slash or backslash
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
  usage-exit "Test group name cannot contain slash/backslash"
else
  test_group=$1
fi

# make directory for today's date, unless it already exists
today=$(date +'%Y-%m-%d')
today_dir="results/$today"

if [ ! -d "$today_dir" ]
then
  echo "================================================================================"
  echo "Creating results directory for "$today" ..."
  mkdir -p "$today_dir"
  echo "Done"
fi

# make main directory for this test group
test_group_dir="$today_dir/$test_group"
  echo "================================================================================"
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

touch $test_group_dir/"middropdata.csv"
touch $test_group_dir/"enddropdata.csv"

if [ ! -f $test_group_dir/$3 ]; then
  echo "Copying config file to test group directory ..."
  cp $3 $test_group_dir
  echo "Done"
fi

# strip comments from config file
file_no_comments=$(sed '/^#/'d $3)

# get param names from header row of config file,
# and save them to an array
param_names=($(echo "$file_no_comments" | awk 'NR==1'))

num_params=${#param_names[@]}

# get number of non-header lines in param file
# sed strips header; xargs trims whitespace
# from return value of wc
num_param_sets=$(echo "$file_no_comments" | sed 1d | wc -l | xargs)

param_set_padding=${#num_param_sets}

# main loop; run the simulation once for each param set
param_set=1
echo "$file_no_comments" | sed 1d | while read -r line; do
  # pad directories to ensure proper ordering
  param_set_dir=$(printf "%s/%0*d" $test_group_dir $param_set_padding $param_set)

  if [ ! -d $param_set_dir ]; then
  echo "================================================================================"
    echo "Creating param set directory: "$param_set_dir" ..."
    mkdir -p $param_set_dir
    echo "Done"
  fi

  # each parameter set will get stored to its own file
  # (as a list of var=val assignments)
  param_set_config_file="$param_set_dir/$test_group-$param_set.conf"

  # save this set of parameter values to an array
  vals=($(echo $line))

  echo "================================================================================"
  echo "Writing parameter set to config file ..."
  i=0
  while [ $i -lt $num_params ]; do
    # construct var=val pair and write it to file
    # Note how crucially this depends on the ordering
    # of param names/values being preserved
    echo "${param_names[i]}=${vals[i]}" >> $param_set_config_file
    let i++
  done

  # echo misc variables to the file
  echo "test_group=$test_group" >> $param_set_config_file
  echo "param_set=$param_set" >> $param_set_config_file
  echo "test_group_dir='$test_group_dir'" >> $param_set_config_file
  echo "param_set_dir='$param_set_dir'" >> $param_set_config_file
  echo "num_param_sets=$num_param_sets" >> $param_set_config_file
  echo "runs_per_param_set=$runs_per_param_set" >> $param_set_config_file
  echo "Done"

  ./run_param_set.sh $param_set_config_file

  let param_set++
done
