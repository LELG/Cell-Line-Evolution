echo "" > middist.dat

for i in $(echo 1 2 3 4 5 6 7 )
do

    cat $1/$i-*/mid_phylo.xml | grep "s:" | grep -v "phylo" | awk '{ print  $2 " " $3 " " $4 }' | tr ":" " " | awk '{print $2 " " $4 " " $6}' | awk '{ if ($3 >=100) print }' | sort -k 2 -n  >> middist.dat

    echo "END" >> middist.dat

done
