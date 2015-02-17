from __future__ import print_function
import argparse
import os
import sys

"""Read output of runs and summarise per group"""

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', default='filename')
parser.add_argument('-l', '--latex', action="store_true", default=False)
parser.add_argument('-p', '--plot', action="store_true", default=False)
opt = parser.parse_args()

total = 0
part = 0
full = 0
maxtime = 0
cmintime = 0
cmaxtime = 0
skip_header = 2
curr_line = "nothing"
avg_precrash_max = 0
avg_postcrash_min = 0
avg_postcrash_max = 0

if opt.latex:
    print("%% full, part, total, mut, sel_pre, prob pos mut, prob neg mut,\
            prob of mut incr, prob of mut decr, min time, max time,\
            cmin time, cmax time")
elif not opt.plot:
    print("FULL, PART, pro,die,mut,sel_pre,prob pos mut, prob neg mut,\
            prob of mut incr, prob of mut decr,\
            maxtime, cmintime, cmaxtime")

lines = [line.strip() for line in open(opt.filename)]
for l in lines:
    if skip_header > 0:
        skip_header -= 1
    else:
        if l == "END_GROUP" and total > 0:
            x = curr_line.split(',')

            if full > 0:
                avg_precrash_max  = maxtime / float(full)
                avg_postcrash_min = cmintime / float(full)
                avg_postcrash_max = cmaxtime / float(full)
            

            if opt.latex:
                print(full, " & ", part, " & ", total, " & ", 
                        x[8],  " & ", #mut
                        x[17][:8],  " & ", #end_mut
                        x[10], " & ",  #sel pressure
                        x[11], " & ",  #prob_mut_pos
                        x[12], " & ",  #prob_mut_neg
                        x[13], " & ",  #prob_inc_mut
                        x[14], " & ",  #prob_dec_mut
                        str(avg_precrash_max)[:6], " & ",  #max_time
                        str(avg_postcrash_min)[:6], " & ",  #cmin_time
                        str(avg_postcrash_max)[:6],  #cmax_time
                        x[32], " & ", #size of all clones that survived crash
                        " \\\\ "
                        )

            elif opt.plot:
                print(float(x[8]), #mut
                      float(x[10]), #sel pre
                      full/float(total))
            else:
                print("FULL ", full, " PART ",part," out of ",total, "details: ",\
                        x[8],x[17][:8],x[10],x[8],x[11],x[12],x[13],x[14],\
                        avg_precrash_max, avg_postcrash_min, avg_postcrash_max, x[32])

            total = 0
            full = 0
            part = 0
            maxtime = 0
            cmaxtime = 0
            cmintime = 0
            avg_precrash_max = 0
            avg_postcrash_min = 0
            avg_postcrash_max = 0
            
        else:
            curr_line = l
            x = l.split(',')
            
            if x[4].strip() == 'FULL' or x[4].strip() == 'FULLNC':
                full += 1
                maxtime += float(x[25])
                cmintime += float(x[27])
                cmaxtime += float(x[29])
                
            if x[4].strip() == 'PART':
                part += 1
            total += 1
