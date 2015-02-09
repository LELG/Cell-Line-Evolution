from __future__ import print_function
import os
import math
from matplotlib import pyplot as plt

def make_popsub(X1, Y1, X2, Y2, filename, title1, title2):

    fig, ax1 = plt.subplots()
    ax1.plot(X1, Y1, color='b', alpha=0.5, linewidth=2)
    ax1.set_xlabel('Discrete Time Intervals', fontsize=9) #color='grey')
    ax1.set_ylabel(title1, fontsize=9, color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
        tl.set_fontsize(9)
    for tl in ax1.get_xticklabels():
        #tl.set_color('grey')
        tl.set_fontsize(9)

    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    ax2 = ax1.twinx()
    ax2.plot(X2, Y2, color='r', alpha=0.5, linewidth=2)
    ax2.set_ylabel(title2, fontsize=9, color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
        tl.set_fontsize(9)

    plt.savefig(filename)
    print("done")
    plt.clf()

def make_plot(X, Y, filename, title):
    #        global DETAILS
    fig, ax1 = plt.subplots()
    ax1.plot(X, Y, color='b', alpha=0.5, linewidth=2)
    ax1.set_xlabel('Discrete Time Intervals', fontsize=9) #color='grey')
    ax1.set_ylabel(title, fontsize=9, color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
        tl.set_fontsize(9)
    for tl in ax1.get_xticklabels():
        #tl.set_color('grey')
        tl.set_fontsize(9)
    #plt.plot( X, Y, '-' )
    #plt.title(title)
    #plt.xlabel( 'Discrete Time Intervals' )
    #plt.ylabel( 'Tumour Size' )
    print("****", filename)
    plt.savefig(filename)
    print("plot done")
    plt.clf()

def make_hist(X, filename, title, bins):
    plt.hist(X, bins, log=True)
    plt.title(title)
    plt.savefig(filename)
    print("plot done")
    plt.clf()

def make_dual_hist(X1, X2, filename, title):
    X1s = []
    X2s = []
    #append for each 'cell'
    for p in X1:
        for i in range(0, p[1]):
            X1s.append(p[0])
    for p in X2:
        for i in range(0, p[1]):
            X2s.append(p[0])

    b = 20

    #plt.hist(X1s, X2s, normed =1, alpha = 0.35, log=True, color=['g','b'])
    plt.hist(X1s, bins=b, alpha=0.45, log=True, linewidth=1)#, normed=1)
    plt.hist(X2s, bins=b, alpha=0.45, log=True, linewidth=1)#, normed=1)
    plt.title(title)
    plt.savefig(filename)
    print("plot done")
    plt.clf()

def make_dual_box(end, mid, filename, title):
    X1s = []
    X2s = []
    #append for each 'cell'
    for p in mid:
        for i in range(0, p[1]):
            X1s.append(p[0])
    for p in end:
        for i in range(0, p[1]):
            X2s.append(p[0])

    data = [X1s, X2s]

    fig, ax1 = plt.subplots()
    bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    plt.title(title)
    plt.savefig(filename)
    plt.clf()

    plt.boxplot(X1s)
    plt.savefig(filename+"boxtest1")
    plt.clf()
    plt.boxplot(X2s)
    plt.savefig(filename+"boxtest2")
    plt.clf()

def make_subpop_life(X, filename, title, end_time, loops, select_time):
    #could 'fatten' or repeat a line if reaches a certain size
    #or colour
    X.sort()
    print("subpop man", filename)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    i = 0
    for col, x1, x2 in X:
        i += 1
        if (x2 == loops+1):
            x2 = end_time
        if (x2 == None):
            x2 = end_time
        if x1 < select_time:
            if col == 'n':
                plt.plot([x1, x2], [i, i], color='green', alpha=0.75)
            else:
                plt.plot([x1, x2], [i, i], color=col)
        else:
            if col == 'n':
                plt.plot([x1, x2], [i, i], color='blue', alpha=0.75)
            else:
                plt.plot([x1, x2], [i, i], color=col)
    ax1.axes.get_yaxis().set_visible(False)
    plt.xlabel('Discrete Time Intervals')
    plt.savefig(filename)
    plt.clf()


def make_subpop_life_mut(X, filename, title, end_time, loops, select_time, mbase, tumoursize):
    #use mutation rate to give colour
    #normalise
    print("subpop man", filename)
    print("base mut rate ", mbase)
    i = 0


    for x1, x2, m, s in X:
        oneprc = tumoursize * 0.01
        while s > 100:
            s = s - 10000
            i += 1
            if (x2 == loops+1):
                x2 = end_time
            if (x2 == None):
                x2 = end_time
            mnorm = m/mbase
            if mnorm < 0:
                mnorm = 0
            if mnorm > 1:
                mnorm = 1
            plt.plot([x1, x2], [i, i], color=(mnorm, 1-mnorm, 0.1))

        i += 1
        if (x2 == loops+1):
            x2 = end_time
        if (x2 == None):
            x2 = end_time
        mnorm = m/mbase
        if mnorm < 0:
            mnorm = 0
        if mnorm > 1:
            mnorm = 1
        plt.plot([x1, x2], [i, i], color=(mnorm, 1-mnorm, 0.1))

    plt.savefig(filename)
    plt.clf()

def mutation_distribution(md, filename, title, SCALE):
    #print("MD",md)
    """ X [mut, precrash size, size] """
    for elem in md:
        if elem[1] > 0:
            plt.scatter(int(elem[1]), int(elem[2]))

    plt.title("Clone size pre and post crash")
    plt.xlabel("Pre-crash size")
    plt.ylabel("Post-crash size")
    plt.savefig(filename)
    plt.clf()

    for elem in md:
        plt.scatter(int(elem[1]), int(elem[2]))

    plt.savefig(filename+"all")
    plt.clf()

def mutation_v_proliferation(clones, filename, title, SCALE):
    plt.axes()
    c2 = 0

    circles = []
    clones.sort(key=lambda tup: tup[2])
    #NEED TO PRE SORT SO LARGER ONES COME OUT ON TOP

    max_clone_size = max(clones,key=lambda x: x[:][2])[2]
    #min_clone_size = min(clones,key=lambda x: x[:][2])[2]

    for c in clones:

        fcol = 'grey'

        c2 = c[2]

        if c2 > 1:
            fcol = 'c'
            if c2 > 50:
                fcol = 'b'
                if c2 > 100:
                    fcol = 'm'
                    if c2 > 200:
                        fcol = 'y'
                        if c2 > 500:
                            fcol = 'r'


        #Get size of clone as ratio of largest
        rad = (c[2] / float(max_clone_size) * 0.0001) + 0.00005 #number could be dynamic
        #rad = 0.0001
        circle = plt.Circle((c[0], c[1]), \
                radius=rad, \
                facecolor='none', edgecolor=fcol)
        plt.gca().add_patch(circle)

    plt.axis('scaled')
    plt.savefig(filename)
    plt.clf()

def mutation_v_proliferation_dat(clones, filename, title, SCALE):

    file_summary = open(filename,'a')

    circles = []
    #NEED TO PRE SORT SO LARGER ONES COME OUT ON TOP
    clones.sort(key=lambda tup: tup[2])
    for c in clones:
        #X - Mut, Y - Pro, Radius - Size
        """
        print(c)
        print(c[0], c[1], c[2],SCALE*10)
        SCALE = SCALE * 10
        c2 = c[2]
        if SCALE == 0:
            SCALE = 10000
        if c[2] == 1:
            c2 = 0.99
        rad =( 1 - (1/math.log(c2)))/float(SCALE) #SCALE*10
        #math.log((c[2])/float(SCALE)), \
        """
        fcol = 'g'

        print(c[0],c[1],c[2],file = file_summary)

def mutation_crash(clones, filename, title, SCALE):
    """ mutation distribution """
    print("mut crash",filename)
    plt.axes()
    c2 = 0

    circles = []
    clones.sort(key=lambda tup: tup[2])
    #NEED TO PRE SORT SO LARGER ONES COME OUT ON TOP

    max_clone_size = max(clones,key=lambda x: x[:][2])[2]
    #min_clone_size = min(clones,key=lambda x: x[:][2])[2]

    for c in clones:

        fcol = 'grey'

        c2 = c[2]

        if c2 > 1:
            fcol = 'c'
            if c2 > 50:
                fcol = 'b'
                if c2 > 100:
                    fcol = 'm'
                    if c2 > 200:
                        fcol = 'y'
                        if c2 > 500:
                            fcol = 'r'
        face = 'none'
        if c[3] > 0:
            fcol = 'g'
            face = 'green'


        #Get size of clone as ratio of largest
        rad = (c[2] / float(max_clone_size) * 0.0001) + 0.00005 #number could be dynamic
        #rad = 0.0001
        circle = plt.Circle((c[0], c[1]), \
                radius=rad, \
                facecolor=face, edgecolor=fcol)
        plt.gca().add_patch(circle)

    plt.axis('scaled')
    plt.savefig(filename)
    plt.clf()
