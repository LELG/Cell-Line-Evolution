from __future__ import print_function
import os
import csv
import subpopulation
import tree_to_xml
import time
import dropdata
from utilities import secs_to_hms

class Analytics():
    """ Record analytics on population

    Population - population size at each point in time
    Subpopulation - number of subpopulations (clones) at each time
    Time - list of time 1... 1000 for plot purposes 
    Mutation - Store average mutation rate at each point in time
    Proliferation - Store average proliferation at each point in time
    Subpopulation Mutation - Store mutation rate for each subpopulation
    Subpopulation Proliferatioin - Store Proliferation rate 
                                    for each subpopulation

        """
    def __init__(self):
        self.population = []
        self.subpopulation = []
        self.time = []
        self.mutation = []
        self.proliferation = []
        self.subpop_mutation = []
        self.subpop_proliferation = []

class Population():
    """ Contains high level simulation routines

        Create initial clone (subpopulation) from input parameters
        Run simulation
        Store data in analytics for each loop
        Print summary data to file
        Print graphs
        
    """
    def __init__(self, opt):
        self.opt = opt
        self.tumoursize = opt.init_size
        self.clonecount = 1
        self.maxsize_lim = opt.maxsize_lim
        self.prolif_lim = self.opt.pro - self.opt.die #could rebase this
        self.opt.prolif_lim = self.opt.pro - self.opt.die
        depth = 0
        time = 0
        mut_type = 'n'
        self.s = subpopulation.Subpopulation(vars(self.opt),\
                 self.opt.pro, self.opt.mut, 
                 depth, time, mut_type, \
                 'n', self.opt.prob_mut_pos, self.opt.prob_mut_neg, self.opt.prob_inc_mut,\
                 self.opt.prob_dec_mut, self.opt.mscale, 0 ) ##I HOPE THIS WORKS
        self.s.size = opt.init_size
        self.analytics_base = Analytics()
        self.select_pressure = 0.0
        self.mutagenic_pressure = 0.0
        self.selective_pressure_applied=False
        for i in range(0,opt.init_diversity):
            self.s.size = 0 #don't use init dummy population if reading from file
            self.s.newsubpop_from_file(self.opt.sub_file)

    def write_population_summary(self, num_cycles, elapsed_time, recovered):
        """ Write details of simulation to master file
        
        In comma delimited format 
        If master file header doesnt exist have bash write it 
        
        """

        # localtime = time.asctime( time.localtime(time.time()) )

        #Get min and max values pre crash 
        min_val, min_time, max_val, max_time = self.precrash_minmax()
        cmin_val = cmin_time = cmax_val = cmax_time = 0

        if self.went_through_crash():
            #Get min and max values post crash 
            #if survived past the crash + buffer
            cmin_val, cmin_time, cmax_val, cmax_time = \
            self.postcrash_minmax()
            #hasty fix for calculting max time 
            if cmax_time == 0 and recovered:
                cmax_time = num_cycles

        recover, recover_type, recover_percent = self.complete_status()

        went_through_crash = 'N'
        size_from_precrash = 0
        if self.went_through_crash():
            went_through_crash = 'Y'
            pp = self.s.tree_to_list("two_side_size")
            #print("pp",pp)
            for (pre,pos) in pp:
                if pre > 0:
                    size_from_precrash+=pos


        #TRACE self.opt.pro

        summary_file = open("{0}/results.dat".format(self.opt.testname), 'a')
        summary_writer = csv.writer(summary_file)

        # assemble values to write
        summary_vals = (self.opt.param_set, self.opt.run_number,
                        went_through_crash,
                        recover, recover_type, recover_percent,
                        self.opt.pro, self.opt.die, self.opt.mut,
                        self.opt.select_time, self.opt.select_pressure,
                        self.s.prob_mut_pos, self.s.prob_mut_neg,
                        self.s.prob_inc_mut, self.s.prob_dec_mut,
                        self.analytics_base.population[-1],    # pop_size
                        self.analytics_base.subpopulation[-1], # num_clones
                        self.analytics_base.mutation[-1],      # avg_mut_rate
                        self.analytics_base.proliferation[-1], # avg_pro_rate
                        secs_to_hms(elapsed_time),
                        len(self.analytics_base.population),   # elapsed_cycles
                        min_val, min_time,
                        max_val, max_time,
                        cmin_val, cmin_time,
                        cmax_val, cmax_time)

        summary_writer.writerow(summary_vals)
        summary_file.close()

    def went_through_crash(self):
        crash_buffer = 25 #check just past crash time
        if len(self.analytics_base.population) > \
        self.opt.select_time + crash_buffer:
            return True
        return False 

    def complete_status(self):
        """ return recovery, full or partial,  percent 

        eg. 'N'  
        eg. 'Y' 


        """
        recover = 'N'
        recover_type = 'NONE'
        
        if self.went_through_crash():
            if self.tumoursize > (self.maxsize_lim/2):
                recover = 'Y'
                recover_type = 'PART'
                if self.tumoursize > (self.maxsize_lim-(self.maxsize_lim/4)): #if > 75% of original size
                    recover_type = 'FULL'
        else: #didnt crash
            if self.tumoursize > (self.maxsize_lim-(self.maxsize_lim/4)):
                recover_type = 'FULLNC'
                recover = 'Y'

        recover_percent = self.tumoursize / float(self.maxsize_lim)

        return recover, recover_type, recover_percent

    def postcrash_minmax(self):
        post_crash_pop =\
                self.analytics_base.population[self.opt.select_time:]
        min_val = min(post_crash_pop) #VALIDATION - can return empty
        min_val_index = post_crash_pop[post_crash_pop.index(min_val)]
        #VALIDATION - can return empty
        #maximum value after lowest point
        max_val = 0
        max_time = 0
        if post_crash_pop[min_val_index:]:
            max_val = max(post_crash_pop[min_val_index:])
            max_time = post_crash_pop.index(max_val) + self.opt.select_time   
        #add back time that we cut out when filtering post crash
        min_time = post_crash_pop.index(min_val) + self.opt.select_time 
        #max time after low point
        return min_val, min_time, max_val, max_time

    def precrash_minmax(self):
        pre_crash_pop =\
                self.analytics_base.population[:self.opt.select_time]
        min_val = min(pre_crash_pop)
        max_val = max(pre_crash_pop)
        min_time = pre_crash_pop.index(min_val) 
        max_time = pre_crash_pop.index(max_val) 
        return min_val, min_time, max_val, max_time

    def info(self):
        print("Parameter Set: ", self.opt)

    def cycle(self,opt):
        self.selective_pressure_applied = False
        extra_lim = self.maxsize_lim * 0.05 #EXTRA LIM NOW 5%
        """ Iteratively cycle through discrete time model

        Recaculate ratio of population size
        Call analytics
        Call subpopulation cycle for each subpopulation
        (Calculate 1 discrete time step for each subpopulation)
        i in loops is TIME

        """

        recovered = False 

        start_time = time.time() # start timing simulation

        for i in range(0, self.opt.loops):
            self.analytics_base.time.append(i)
            ratio = float(self.tumoursize) / float(self.maxsize_lim)
            #beta(alpha=3,beta=1)
            """ NEW BEGIN """
            """
            
            a = 3
            b = 1
            ratio = (ratio**(a-1) * (1-ratio)**(b-1))
            #scale growth non-linearly
           
            """
            """ NEW END """
            prolif_adj = ratio * float(self.prolif_lim)

            #self.s.init_prolif(prolif_adj) #PUSHED TO CYCLE

#CYCLE SHOULD RETURN POP COUNT
            #HAVE cycle return tumoursize
            self.tumoursize, self.clonecount, avg_mut, avg_pro \
                                             = self.s.cycle(self.tumoursize, \
                                               self.select_pressure, 
                                               self.mutagenic_pressure, i, 
                                               prolif_adj) 
            #return 'changes' / number of cells instead of counting each time
            avg_mut_rate = 0
            avg_pro_rate = 0
            if self.tumoursize > 0:
                avg_mut_rate = avg_mut / float(self.tumoursize)
                avg_pro_rate = avg_pro / float(self.tumoursize)
            self.update_analytics(avg_mut_rate, avg_pro_rate)
            if self.opt.A: #auto dynamic restriction of size
                if self.analytics_base.proliferation[-1] > self.opt.pro:
                    #if avg pro higher than normal, use it as top
                    self.prolif_lim = self.analytics_base.proliferation[-1] - self.opt.die
                else:
                    self.prolif_lim = self.opt.pro
            
            if not self.selective_pressure_applied:
                if i == self.opt.select_time:
                    self.s.set_precrash_size()
                    self.selective_pressure()
                    if not self.opt.NP:
                        self.print_results("mid",i)
                    tree_to_xml.tree_parse(self.s, self.tumoursize, i, "mid0")
                    if opt.init_diversity:
                        dropdata.drop(self.s, self.tumoursize, i, "mid0")
                    self.selective_pressure_applied = True
                    #PRINT RESULTS with diff filename

                #If maxsize flag is set auto start selective pressure
                #at size lim
                else:
                    if self.opt.M:
                        if self.tumoursize > self.maxsize_lim:
                            self.s.set_precrash_size()
                            self.selective_pressure()
                            if not self.opt.NP:
                                self.print_results("mid",i)
                            tree_to_xml.tree_parse(self.s, self.tumoursize, i, "mid")
                            if opt.init_diversity:
                                dropdata.drop(self.s, self.tumoursize, i, "mid")
                            self.opt.select_time = i #update time when sel press introduced
                            self.selective_pressure_applied = True
                            #PRINT RESULTS with diff filename


            #PRINT INFO EVERY 100 LOOPS?
            if i % 1000 == 0:
                os.system('clear')
                print(i, " out of", self.opt.loops, " tumour size: ", \
                    self.tumoursize)
            if self.tumoursize == 0:
                if not self.selective_pressure_applied:
                    return 0
                break
            if self.tumoursize > (self.maxsize_lim + extra_lim):
                #IF TOO BIG BUT HASNT CRASHED, LET NEXT CYCLE CRASH IT
                if self.selective_pressure_applied:
                    print("TUMOUR IS TOO BIG")
                    recovered = True
                    break

            #self.s.info()

        print("cycle done")
        end_time = time.time() # finish timing simulation
        elapsed_time = end_time - start_time
        self.write_population_summary(i, elapsed_time, recovered)

        if not self.opt.NP:
            self.print_results("end",i)
            self.print_plots("new",i)

        #ALWAYS PLOT TREE
        fname=""
        tree_to_xml.tree_parse(self.s,self.tumoursize,i,fname)
        if opt.init_diversity:
            dropdata.drop(self.s, self.tumoursize, i, "end")
            print("printing the DROP")
        return 1
        #print(self.analytics_base.population)

    def selective_pressure(self):
        self.select_pressure = self.opt.select_pressure
        self.mutagenic_pressure = self.opt.mutagenic_pressure
        self.mid_proliferation = self.s.tree_to_list("proliferation_size")
        self.mid_mutation = self.s.tree_to_list("mutation_rate")

    def update_analytics(self,avg_mut_rate,avg_pro_rate):
        #print("analytics")
        sub = self.clonecount
        self.analytics_base.population.append(self.tumoursize)
        self.analytics_base.subpopulation.append(sub)

        #EFFECTIVE PROLIFERATION
        if self.analytics_base.time[-1] > self.opt.select_time:
            self.analytics_base.proliferation.append(\
                    avg_pro_rate - self.opt.select_pressure)
        else:
            self.analytics_base.proliferation.append(\
                    avg_pro_rate) 

        self.analytics_base.mutation.append(avg_mut_rate)

    def print_results(self,when,end_time):
        """ Print all results to plots / file


        Print result to graphs using matplotlib
        If --R option parsed print raw output to file

        """

        filename = self.opt.filename+when

        if self.opt.R:
            from outputdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat
        else:
            from plotdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat,\
                    mutation_distribution,  mutation_crash, make_subpop_life_mut,\
                    make_popsub

        """ Print only once at end of simulation """

        if self.selective_pressure_applied: #ONLY PRINT AT END OF SIM
            #POPULATION
            make_plot(self.analytics_base.time,
                    self.analytics_base.population,
                    filename+"population_graph","Population Size") 
            #SUBPOPULATION
            make_plot(self.analytics_base.time,
                    self.analytics_base.subpopulation,
                    filename+"subpop_graph","No. of Clones") 

            #POPULATION + SUBPOPULATION
            make_popsub(self.analytics_base.time, self.analytics_base.population,
                    self.analytics_base.time, self.analytics_base.subpopulation,
                    filename+"popsubpop",'Tumour Size','No. of Clones')
            #MUTATION RATE BUT SUBPOPULATION
            make_popsub(self.analytics_base.time, self.analytics_base.mutation,
                    self.analytics_base.time, self.analytics_base.subpopulation,
                    filename+"mutsubpop",'Average Mutation Rate','No. of Clones')

            #PROLIFERATION RATE + MUTATION RATE
            make_popsub(self.analytics_base.time, self.analytics_base.mutation,
                    self.analytics_base.time, self.analytics_base.proliferation,
                    filename+"prolifmut",'Mutation Rate','Proliferation Rate')
            #PROLIFERATION RATE + POPULATION
            make_popsub(self.analytics_base.time, self.analytics_base.population,
                    self.analytics_base.time, self.analytics_base.proliferation,
                    filename+"prolifandpop",'Tumour Size','Proliferation Rate')

            #EFFECTIVE PROLIFERATION RATE
            make_plot(self.analytics_base.time,
                    self.analytics_base.proliferation,
                    filename+"effect_prolif","EFFECT PROLIFERATION")
            #MUTATION RATE AVG
            make_plot(self.analytics_base.time, self.analytics_base.mutation,
                  filename+"mutation_avg","MUTATION AVG")

            # Mutation...
            md = self.s.tree_to_list("mutation_distribution")
            mutation_distribution(md, \
                    filename+"mutation_distribution", \
                    "MUTATION V TIME - PRE/POST CRASH",\
                    self.opt.scale)

            # Mutation
            if (self.tumoursize > 0):
                md = self.s.tree_to_list("mutation_distribution_1")
                mutation_crash(md, \
                        filename+"mutation_distribution_1", \
                        "MUTATION V TIME - PRE/POST CRASH",\
                        self.opt.scale)

            #cell lines graph  - [(self.s_time,self.d_time)]
            make_subpop_life(self.s.tree_to_list("cell_line_time"),
                    filename+"cell_lines_alpha",
                    "CELL LIFESPAN",end_time,self.opt.loops, 
                    self.opt.select_time)

            """
            #cell lines graph  - [(self.s_time,self.d_time)]
            make_subpop_life_mut(self.s.tree_to_list("cell_line_time_mut"),
                    filename+"cell_lines_muta",
                    "CELL LIFESPAN",end_time,self.opt.loops, 
                    self.opt.select_time,self.opt.mut,self.tumoursize)
            """


        """ Print at mid and end of simulation """       

        #PRINT END HISTOGRAM IF POPULATION STILL ALIVE
        if (self.tumoursize > 0):
            #PROLIFERATION HISTOGRAM
            # [(self.proliferation-self.prolif_adj,self.size)]

            pro_hist = self.s.tree_to_list("proliferation")
            pro_hist.sort()
            mut_hist = self.s.tree_to_list("mutation")
            mut_hist.sort()

            make_hist(pro_hist,
                      filename+"proliferation_hist",
                      "PROLIFERATION RATES",25) 
            #MUTATION HISTOGRAM
            # [self.mutation]
            make_hist(mut_hist,
                      filename+"mutation_hist",
                      "MUTATION RATES",25) 
            #POPULATION HISTOGRAM
            # [self.size]

            pop_hist = self.s.tree_to_list("size"),

            make_hist(pop_hist,
                      filename+"population_hist",
                      "POPULATION DIVISION",25)
            #ALLELE FREQ
            # just_allele_freq z = z + [i/float(tumoursize)]
            norm, r, just_allele_freq = \
                    self.s.freq_of_mutation(self.tumoursize)
            make_hist(just_allele_freq,
                    filename+"allele", "ALLELE FREQ",self.analytics_base.subpopulation[-1]) #bins equal to number of sub pops
            #CELL CIRCLE MUT V PRO RATES
            # if size > 0 [(self.mutation,self.proliferation,self.size)]
            mutation_v_proliferation(self.s.tree_to_list("circles"), \
                    filename+"circles", \
                    "MUTATION V PROLIFERATION RATES",\
                    self.opt.scale)

            # [(self.mutation,self.proliferation,self.size)]
            mutation_v_proliferation(self.s.tree_to_list("circles_all"), \
                    filename+"circles_all", \
                    "MUTATION V PROLIFERATION RATES",\
                    self.opt.scale)

            #MAKE CIRCLES ACROSS ALL GRAPHS BY WRITING TO 1 FILE
            # if size > 0 [(self.mutation,self.proliferation,self.size)]
            mutation_v_proliferation_dat(self.s.tree_to_list("circles"), \
                    self.opt.testgroup+"circles.dat", \
                    "MUTATION V PROLIFERATION RATES",\
                    self.opt.scale)
            # [(self.mutation,self.proliferation,self.size)]
            mutation_v_proliferation_dat(self.s.tree_to_list("circles_all"), \
                    self.opt.testgroup+"circles_all.dat", \
                    "MUTATION V PROLIFERATION RATES",\
                    self.opt.scale)


    def print_plots(self,when,end_time):
        """ Print all results to plots / file


        Print result to graphs using matplotlib
        If --R option parsed print raw output to file

        """

        filename = self.opt.filename+when

        if self.opt.R:
            from outputdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat
        else:
            from plotdata import make_plot, make_subpop_life, make_hist, \
                    mutation_v_proliferation, mutation_v_proliferation_dat,\
                    mutation_distribution,  mutation_crash, make_subpop_life_mut, \
                    make_dual_hist, make_dual_box

        """ Proliferation Histogram """

        if (self.tumoursize > 0):

            end_proliferation = self.s.tree_to_list("proliferation_size")
            end_mutation = self.s.tree_to_list("mutation_rate")

            #PROLIFERATION HISTOGRAM
            # [(self.proliferation-self.prolif_adj,self.size)]

            end_proliferation.sort()
            self.mid_proliferation.sort()
            make_dual_hist(end_proliferation, self.mid_proliferation,
                        filename+"prolif_hist", "Proliferation Rates")
            
            end_mutation.sort()
            self.mid_mutation.sort()
            make_dual_hist(end_mutation, self.mid_mutation,
                        filename+"mutation_hist", "Mutation Rates")
            
            make_dual_box(end_mutation, self.mid_mutation,
                        filename+"mutation_box", "Mutation Rate Box")
