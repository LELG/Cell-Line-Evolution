# filter those results > 100, and show how many left
# then if 100 tests, sum of results / 100,000 *100= ~99.4% of results this accounts for
cat ZERO_INC/*/_phylo.xml | grep "s:" | grep -v "phylo" | awk '{ print $3 " " $4 " " $5 }' | tr ":" " " |     awk '{print $2 " "$4 " " $6}' | grep -v ' n' | sort -k 2 -n | awk '{ if ($2 >= 100) print }' | awk '{sum[$3] += $2;}END{for (s in sum){print sum[s], s;}}'
