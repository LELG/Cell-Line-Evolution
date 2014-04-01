cat $1/*/_phylo.xml | grep "s:" | grep -v "phylo" | awk '{ print $2 " " $3 " " $4 " " $5 }' | tr ":" " " | awk '{print $2 " "$4 " " $6}' | sort -k 3 -n | awk '{ if ($3 >=100) print }' > dist.dat
#python distribution_plot.py -f dist.dat -o plots/$2
