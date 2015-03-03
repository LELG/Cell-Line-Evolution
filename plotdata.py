"""
Functions for plotting simulation data.

Authors
-------
Andrew Bakshi : andrew.bakshi@gmail.com
Yoshua Wakeham : yoshwakeham@gmail.com
"""

from __future__ import print_function
from matplotlib import pyplot as plt

def make_dual_plot(xdata, y1data, y2data, filename, title1, title2):
    """Plot two dependent vars against the same independent var."""
    _fig, ax1 = plt.subplots()
    ax1.plot(xdata, y1data, color='b', alpha=0.5, linewidth=2)
    ax1.set_xlabel('Discrete Time Intervals', fontsize=9)
    ax1.set_ylabel(title1, fontsize=9, color='b')
    for label in ax1.get_yticklabels():
        label.set_color('b')
        label.set_fontsize(9)
    for label in ax1.get_xticklabels():
        label.set_fontsize(9)

    ax2 = ax1.twinx()
    ax2.plot(xdata, y2data, color='r', alpha=0.5, linewidth=2)
    ax2.set_ylabel(title2, fontsize=9, color='r')
    for label in ax2.get_yticklabels():
        label.set_color('r')
        label.set_fontsize(9)

    plt.savefig(filename)
    print("PLOT CREATED: " + filename)
    plt.close()

def make_plot(xdata, ydata, filename, title):
    """Plot a single independent variable."""
    #        global DETAILS
    _fig, ax1 = plt.subplots()
    ax1.plot(xdata, ydata, alpha=0.5, linewidth=2)
    ax1.set_xlabel('Discrete Time Intervals', fontsize=9)
    ax1.set_ylabel(title, fontsize=9, color='b')
    for label in ax1.get_yticklabels():
        label.set_color('b')
        label.set_fontsize(9)
    for label in ax1.get_xticklabels():
        label.set_color('grey')
        label.set_fontsize(9)
    plt.savefig(filename)
    print("PLOT CREATED: " + filename)
    plt.close()

def make_hist(xdata, filename, title, bins=25, log=False):
    """Plot a histogram."""
    plt.hist(xdata, bins, log=log)
    plt.title(title)
    plt.savefig(filename)
    print("HISTOGRAM CREATED: " + filename)
    plt.close()

def make_dual_hist(x1freqs, x2freqs, filename, title):
    """Plot two histograms on the same axes."""
    x1data = []
    x2data = []

    #append for each 'cell'
    for x1var, freq in x1freqs:
        for _ in range(freq):
            x1data.append(x1var)
    for x2var, freq in x2freqs:
        for _ in range(freq):
            x2data.append(x2var)

    num_bins = 20

    plt.hist(x1data, bins=num_bins,
             alpha=0.45, log=True, linewidth=1)
    plt.hist(x2data, bins=num_bins,
             alpha=0.45, log=True, linewidth=1)
    plt.title(title)
    plt.savefig(filename)
    print("HISTOGRAM CREATED: " + filename)
    plt.close()

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
    plt.close()

def make_subpop_life(X, filename, title, end_time, loops, select_time):
    #could 'fatten' or repeat a line if reaches a certain size
    #or colour
    X.sort()
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
    print("PLOT CREATED: ", filename)
    plt.close()


def make_subpop_life_mut(X, filename, title, end_time, loops, select_time, mbase, tumoursize):
    #use mutation rate to give colour
    #normalise
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
    print("PLOT CREATED: ", filename)
    print("base mut rate ", mbase)
    plt.close()

def plot_mut_effect_sizes(mutations, rate_type, filename, bins=100):
    """Plot the distribution of mutation effect sizes."""
    if rate_type == 'prolif':
        effect_data = [m.prolif_rate_effect for m in mutations]
    elif rate_type == 'mut':
        effect_data = [m.mut_rate_effect for m in mutations]
    else:
        raise ValueError("Rate type must be `prolif` or `mut`")

    title = "Mutation Effect Sizes ({} Rate)".format(rate_type.title())
    make_hist(effect_data, filename, title, bins)

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
    print("PLOT CREATED: " + filename + "all")
    plt.close()

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
    print("PLOT CREATED: ", filename)
    plt.close()

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
    print("PLOT CREATED: ", filename)
    plt.close()


def print_results(popn, when, end_time, plot_style=None):
    """ Print all results to plots / file

    Print result to graphs using matplotlib
    If --r_output option parsed print raw output to file

    """
    if plot_style:
        plt.style.use(plot_style)

    filename = "{0}/{1}_".format(popn.opt.run_dir, when)

    anlt = popn.analytics_base

    """
    if popn.opt.r_output:
        from outputdata import make_plot, make_subpop_life, make_hist, \
                mutation_v_proliferation, mutation_v_proliferation_dat
    """

    # Only make these plots at the end of the sim
    if popn.selective_pressure_applied:
        # Population vs Time
        make_plot(anlt.time, anlt.tumoursize,
                  filename + "population_graph", "Population Size")
        # Clone count vs Time
        make_plot(anlt.time, anlt.clonecount,
                  filename + "subpop_graph", "No. of Clones")
        # Effective Proliferation Rate
        make_plot(anlt.time, anlt.avg_proliferation,
                  filename + "effect_prolif", "Effective Proliferation Rate")
        # Avg Mutation Rate
        make_plot(anlt.time, anlt.avg_mutation,
                  filename + "mutation_avg", "Average Mutation Rate")
        # Selective pressure
        make_plot(anlt.time, anlt.select_pressure,
                  filename + "select_pressure", "Qty of Selective Pressure")

        # Population vs Clone count
        make_dual_plot(anlt.time,
                       anlt.tumoursize, anlt.clonecount,
                       filename + "popsubpop",
                       'Tumour Size', 'No. of Clones')
        # Mutation rate vs Clone count
        make_dual_plot(anlt.time,
                       anlt.avg_mutation, anlt.clonecount,
                       filename + "mutsubpop",
                       'Average Mutation Rate', 'No. of Clones')
        # Proliferation Rate vs Mutation Rate
        make_dual_plot(anlt.time,
                       anlt.avg_mutation, anlt.avg_proliferation,
                       filename + "prolifmut",
                       'Mutation Rate', 'Proliferation Rate')
        # Proliferation Rate vs Population
        make_dual_plot(anlt.time,
                       anlt.tumoursize, anlt.avg_proliferation,
                       filename + "prolifandpop",
                       'Tumour Size', 'Proliferation Rate')
        # Selective pressure vs Population
        make_dual_plot(anlt.time,
                       anlt.tumoursize, anlt.select_pressure,
                       filename + "pop_v_select_pressure",
                       'Tumour Size', 'Selective Pressure')

        # Histogram of mutation effect sizes
        plot_mut_effect_sizes(popn.all_mutations, 'prolif',
                              filename + "mut_effect_sizes_prolif")
        plot_mut_effect_sizes(popn.all_mutations, 'mut',
                              filename + "mut_effect_sizes_mut")

        # Mutation...
        mut_distro = popn.subpop.tree_to_list("mutation_distribution")
        mutation_distribution(mut_distro,
                              filename + "mutation_distribution",
                              "Mutation vs Time - Pre/Post Crash",
                              popn.opt.scale)

        # Mutation
        if not popn.is_dead():
            mut_distro = popn.subpop.tree_to_list("mutation_distribution_1")
            mutation_crash(mut_distro,
                           filename + "mutation_distribution_1",
                           "Mutation vs Time - Pre/Post Crash",
                           popn.opt.scale)

        #cell lines graph  - [(popn.s_time,popn.d_time)]
        make_subpop_life(popn.subpop.tree_to_list("cell_line_time"),
                         filename + "cell_lines_alpha",
                         "Cell Lifespan", end_time, popn.opt.max_cycles,
                         popn.opt.select_time)

    # Print these plots both at the crash, and at the end of sim
    if not popn.is_dead():
        #PROLIFERATION HISTOGRAM
        # [(popn.proliferation-popn.prolif_adj,popn.size)]

        pro_hist = popn.subpop.tree_to_list("proliferation")
        pro_hist.sort()
        mut_hist = popn.subpop.tree_to_list("mutation")
        mut_hist.sort()

        make_hist(pro_hist,
                  filename+"proliferation_hist",
                  "Proliferation Rates", log=True)
        #MUTATION HISTOGRAM
        # [popn.mutation]
        make_hist(mut_hist,
                  filename+"mutation_hist",
                  "Mutation Rates", log=True)
        #POPULATION HISTOGRAM
        # [popn.size]

        pop_hist = popn.subpop.tree_to_list("size"),

        make_hist(pop_hist,
                  filename+"population_hist",
                  "Population Division", log=True)
        #ALLELE FREQ
        # just_allele_freq z = z + [i/float(tumoursize)]
        norm, r, just_allele_freq = \
                popn.subpop.freq_of_mutation(popn.tumoursize)
        make_hist(just_allele_freq,
                  filename + "allele",
                  "Allele Freq",
                  #bins equal to number of sub pops
                  bins=anlt.clonecount[-1], log=True)
        #CELL CIRCLE MUT V PRO RATES
        # if size > 0 [(popn.mutation,popn.proliferation,popn.size)]
        mutation_v_proliferation(popn.subpop.tree_to_list("circles"),
                                 filename + "circles",
                                 "Mutation vs Proliferation Rates",
                                 popn.opt.scale)

        # [(popn.mutation,popn.proliferation,popn.size)]
        mutation_v_proliferation(popn.subpop.tree_to_list("circles_all"),
                                 filename + "circles_all",
                                 "Mutation vs Proliferation Rates",
                                 popn.opt.scale)

        #MAKE CIRCLES ACROSS ALL GRAPHS BY WRITING TO 1 FILE
        # if size > 0 [(popn.mutation,popn.proliferation,popn.size)]
        fpath = "{0}/{1}-circles.dat".format(popn.opt.param_set_dir,
                                             popn.opt.param_set)
        mutation_v_proliferation_dat(popn.subpop.tree_to_list("circles"),
                                     fpath,
                                     "Mutation vs Proliferation Rates",
                                     popn.opt.scale)
        # [(popn.mutation,popn.proliferation,popn.size)]
        fpath = "{0}/{1}-circles_all.dat".format(popn.opt.param_set_dir,
                                                 popn.opt.param_set)
        mutation_v_proliferation_dat(popn.subpop.tree_to_list("circles_all"),
                                     fpath,
                                     "Mutation vs Proliferation Rates",
                                     popn.opt.scale)


def print_plots(popn, when):
    """ Print all results to plots / file

    Print result to graphs using matplotlib
    If --r_output option parsed print raw output to file
    """

    filename = "{0}/{1}_".format(popn.opt.run_dir, when)

    """
    if popn.opt.r_output:
        from outputdata import make_plot, make_subpop_life, make_hist, \
                mutation_v_proliferation, mutation_v_proliferation_dat
    """

    # Proliferation Histogram #

    if not popn.is_dead():
        end_proliferation = popn.subpop.tree_to_list("proliferation_size")
        end_mutation = popn.subpop.tree_to_list("mutation_rate")

        #PROLIFERATION HISTOGRAM
        # [(popn.proliferation-popn.prolif_adj,popn.size)]

        end_proliferation.sort()
        popn.mid_proliferation.sort()
        make_dual_hist(end_proliferation, popn.mid_proliferation,
                       filename + "prolif_hist", "Proliferation Rates")

        end_mutation.sort()
        popn.mid_mutation.sort()
        make_dual_hist(end_mutation, popn.mid_mutation,
                       filename + "mutation_hist", "Mutation Rates")

        make_dual_box(end_mutation, popn.mid_mutation,
                      filename + "mutation_box", "Mutation Rate Box")
