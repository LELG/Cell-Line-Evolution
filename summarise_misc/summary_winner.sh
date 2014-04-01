cat $1/results.dat | grep FULL | awk '{ print $3 }' > mut_recovered.dat 

echo "" > blank_win.txt
rm blank_win.txt

for i in $(cat mut_recovered.dat)
do
    cat $i/_phylo.xml | grep "s:" | grep -v "phylo" | awk '{ print $4 " " $5 }' | tr ":" " " | awk '{ print $2 " " $4 }' | awk '{ print $1 " " $2 }' | awk '{sum[$2] += $1;}END{for (s in sum){print sum[s], s;}}' | sort -k 1 -n | tail -1 >> blank_win.txt
    #echo "survived till end: "
    #cat $i/_phylo.xml | grep "s:" | grep -v "phylo" | awk '{ print $4 " " $5 }' | tr ":" " " | awk '{ print $2 " " $4 }' | awk '{ print $1 " " $2 }' | awk '{sum[$2] += $1;}END{    for (s in sum){print sum[s], s;}}' | grep -v "0 n"
    

    #echo "SUMMARY:"

done

python distribution_wins.py -f blank_win.txt -o plots/$2
