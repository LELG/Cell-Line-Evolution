#just 100 runs of low baseline mutation rates

if [ -z $1 ]
then
	echo "USAGE: "$0" [directory name] [test name]"
	echo "Need directory as first command line argument"
	exit
fi

dirname=$1
testname=$1


if [ ! -d $dirname ]
then
	mkdir $dirname
fi

currentdir=$(pwd)

## NP - no plots (graphs)

proliferation_rate=0.04
death_rate=0.03
mutation_rate=0.001
loops=1000000
maxsize_lim=100000
initial_size=25
initial_diversity=0
selective_pressure=0.01
select_time=400000
prob_mut_pos=0.01
prob_mut_neg=0.99 #zero neutral
prob_inc_mut=0.0
prob_dec_mut=0.0
test_per_permuation=100 #run each 0 to n times
scale=0.5
mscale=10.0
sub_file="mut.dat"

#FLAGS: --M auto selective pressure @ max size

#mutation_values="0.0000001 0.0000005 0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05"
mutation_values="0.0000001 0.000001 0.00001 0.0001 0.001 0.005 0.01"
incr_mut_rate="0.99"

pre_header="mutation_values: "$mutation_values", tests per group: "$test_per_permutation

header=">,filename,went through crash,recovered,recovery type, recovery percent, pro, die, mut, select time, select pressyre, prob pos mut, prob neg mut, prob of increase in mut rate, prob of mut incr, prob of mut decr, pop size, no. of clones, avg mut rate at end, avg pro rate at end, time, tumoursize, cycles, pre crash in, pre crash min time, pre crash max, ax time, post crash min, min time, post crash max, max time"

#mkdir -p $testname/results/

echo $pre_header >> $testname/"results.dat"
echo $header >> $testname/"results.dat"
echo "TESTSIZE, etc"

count=1

#MULTIPLE TESTS
for mutation_rate in $mutation_values
do
    for prob_inc_mut in $incr_mut_rate
    do

        prob_dec_mut=$(echo $prob_inc_mut | awk '{ print 1-$1}')

        j=0
        while [ $j -lt $test_per_permuation ]
        do

            j=$((j+1))

            filepath=$dirname/$count-$j/
            testgroup=$dirname/$count

            if [ ! -d $filepath ]
            then
                mkdir $filepath
            fi

            python main.py -d $death_rate -p $proliferation_rate -m $mutation_rate --loops $loops --maxsize_lim $maxsize_lim --prolif_lim 0.0 --init_size $initial_size -f $filepath -s $selective_pressure -t $select_time --prob_mut_pos $prob_mut_pos --prob_mut_neg $prob_mut_neg --prob_inc_mut $prob_inc_mut --prob_dec_mut $prob_dec_mut --scale $scale --mscale $mscale --M -n $testname -g $testgroup --init_diversity $initial_diversity --sub_file $sub_file --Z --NP

        done
        count=$(($count+1))

    #INDICATE AN IDENTICAL SET OF TESTS FINISHED
    echo "END_GROUP" >> $testname/"results.dat"
    #python plot_circles.py -f $testgroup"_circles.dat"
    #python plot_circles.py -f $testgroup"_circles_all.dat"
    done

done


#python compact.py -f $testname/results.dat > $testname/summary.dat

echo $testgroup
