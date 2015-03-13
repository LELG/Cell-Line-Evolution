"""
Module for classes and functions related to
individual clones.

Authors
-------
Andrew Bakshi : andrew.bakshi@gmail.com
Yoshua Wakeham : yoshwakeham@gmail.com
"""
from __future__ import print_function
import random
import csv
import json
import numpy as np
from mutation import Mutation

class Subpopulation(object):
    """
    Individual cancer clone.
    """
    def __init__(self, opt, prolif, mut_rate, depth, t_curr, col, prev_time, inherited_mutns=None):
        """Create new clone."""
        self.prolif_rate = prolif
        self.mut_rate = mut_rate
        # self.mut_static = m #static for scale calution of new value
        self.death_rate = opt.die
        self.size = 1
        self.precrash_size = 0
        self.nodes = []
        self.depth = depth
        self.s_time = t_curr
        self.d_time = None
        #self.idnt = "" ##This gets created when we run freq through tree
        self.col = col
        self.branch_length = t_curr - prev_time

        if inherited_mutns is None:
            self.mutations = {'b': [], 'n': [], 'd': [], 'r': []}
        else:
            self.mutations = inherited_mutns

        if len(self.mutations['r']) > 0:
            self.is_resistant = True
        else:
            self.is_resistant = False

    def new_subpop_from_file(self, opt, filename):
        """Load a heterogeneous starting population from a subfile."""
        subfile = open(filename)
        reader = csv.DictReader(subfile)
        for line in reader:
            mut = float(line["#mut"])
            col = line["col"]

            #add new subpopulation node
            new_subpop = Subpopulation(opt=opt,
                                       prolif=self.prolif_rate,
                                       mut_rate=mut, depth=1, t_curr=0,
                                       col=col, prev_time=0)
            new_subpop.size = opt.init_size
            self.nodes.append(new_subpop)

    def info(self):
        """Print info about this clone and all its children"""
        print("Clone Info")
        print(self.__dict__)
        print("No. of nodes: ", len(self.nodes))
        print(self.prolif_rate, self.mut_rate)
        for node in self.nodes:
            node.info()

    def set_precrash_size(self):
        """ Record pre crash size """
        self.precrash_size = self.size
        for node in self.nodes:
            node.set_precrash_size()

    def update(self, opt, select_pressure, mutagenic_pressure, t_curr, prolif_adj, all_muts):
        """Update this clone and its children for one time step."""
        if self.is_dead_end():
            return (0, 0, 0, 0,)

        new_pop_size = new_sub_count = new_mut_agg = new_pro_agg = 0

        if not self.is_dead():
            effective_prolif = self.prolif_rate - prolif_adj - select_pressure
            if mutagenic_pressure:
                effective_mut = self.mut_rate * mutagenic_pressure
            else:
                effective_mut = self.mut_rate
            # sample for cell death, division and mutation
            # note that we sample for both division AND death before
            # updating the clone size. This means that a 'cell'
            # could technically reproduce and die in the same cycle
            # (i.e. if cells_dead + cells_new > initial_size)
            cells_new = safe_binomial_sample(self.size, effective_prolif)
            cells_dead = safe_binomial_sample(self.size, self.death_rate)
            self.size = self.size + cells_new - cells_dead
            # this is the total number of mutations this cycle,
            # not necessarily number of new subclones to spawn
            new_mutns = safe_binomial_sample(cells_new, effective_mut)
        else:
            # clone is dead - update its attributes
            # if it died on the previous cycle
            if self.size < 0:
                self.size = 0
            if not self.d_time:
                self.d_time = t_curr

        # update child nodes - whether or not clone is alive
        for node in self.nodes:
            node_results = node.update(opt, select_pressure, mutagenic_pressure,
                                       t_curr, prolif_adj, all_muts)
            node_pop, node_sub_count, node_mut_agg, node_pro_agg = node_results
            new_pop_size += node_pop
            new_sub_count += node_sub_count
            new_mut_agg += node_mut_agg
            new_pro_agg += node_pro_agg

        # check again whether clone is alive;
        # clones which have died in this update
        # will now register as dead
        if not self.is_dead():
            for _i in xrange(new_mutns):
                new_mutn = Mutation(opt, all_muts)
                self.new_child(t_curr, opt, new_mutn)
                self.size -= 1
                new_sub_count += 1
                new_pop_size += 1

        # finally, add this clone's stats to the return values
        new_pop_size += self.size
        new_sub_count += 1
        # return aggregate prolif and mut rates, ignoring
        # selective and mutagenic pressure
        new_mut_agg += self.mut_rate * self.size
        new_pro_agg += (self.prolif_rate - prolif_adj) * self.size

        return new_pop_size, new_sub_count, new_mut_agg, new_pro_agg

    def new_child(self, t_curr, opt, new_mutn):
        """Spawn a new child clone."""
        max_pro_rate = 1.0
        max_mut_rate = 1.0
        min_pro_rate = 0.0
        min_mut_rate = 1e-10

        # get new prolif rate, making sure we
        # bound above and below
        provis_prolif = self.prolif_rate + new_mutn.prolif_rate_effect
        new_prolif_rate = max(min_pro_rate, min(max_pro_rate, provis_prolif))

        # get new mutation rate, again, bounding appropriately
        provis_mut = self.mut_rate + new_mutn.mut_rate_effect
        new_mut_rate = max(min_mut_rate, min(max_mut_rate, provis_mut))

        # create shallow copy of mutation dictionary
        child_mutns = self.mutations.copy()
        # now create shallow copies of each mut type list.
        # the lists will still point to identical mutn objects,
        # but appending to parent's list will not change child's list,
        # and vice versa
        for mut_type in child_mutns:
            child_mutns[mut_type] = self.mutations[mut_type][:]

        # now append new mutation to child's inherited mutations
        child_mutns[new_mutn.mut_type].append(new_mutn)

        new_depth = self.depth + 1

        child = Subpopulation(opt=opt,
                              prolif=new_prolif_rate, mut_rate=new_mut_rate,
                              depth=new_depth, t_curr=t_curr,
                              col=self.col, prev_time=self.s_time,
                              inherited_mutns=child_mutns)

        new_mutn.original_clone = child
        self.nodes.append(child)

    def switch_mutn_type(self, mutn, new_type):
        """Change the type of one of this clone's mutations."""
        try:
            self.mutations[mutn.mut_type].remove(mutn)
        except ValueError:
            raise
        try:
            self.mutations[new_type].append(mutn)
        except KeyError:
            raise

        for node in self.nodes:
            node.switch_mutn_type(mutn, new_type)

    def become_resistant(self):
        """
        Register the fact that this clone contains a resistance mutation.
        """
        self.is_resistant = True
        for node in self.nodes:
            node.become_resistant()

    def prune_dead_end_clones(self):
        """Remove all dead end clones from subpopulation tree.

        A dead end clone is defined as one
        which is dead, and has no children.
        """
        children_to_keep = [node for node
                            in self.nodes
                            if not node.is_dead_end()]
        self.nodes = children_to_keep

        for node in self.nodes:
            node.prune_dead_end_clones()


    # GATHER DATA FOR ANALYTICS
    # Dangerous loops beyond this sign

    def pop_as_list(self):
        n = self.get_clone_attrs_as_list("size")
        n.sort()
        #out = filter(lambda a: a!=0, n)
        return n

    def pop_count(self):
        x = self.pop_as_list()
        z = 0 
        for i in x:
            z = z + i
        return z

    def freq_of_mutation(self,tumoursize):
        x = self.freq_to_list("0-0")
        y = []
        z = []
        for i, j in x:
            y = y + [i/float(tumoursize), j]
            z = z + [i/float(tumoursize)]
        return x, y, z

    def freq_to_list(self, idnt):
        self_node = []
        blank_ctree = []
        # TODO fix this - it's currently broken.
        # TODO Should subpop store type of its initial mutation?
        #if self.mutations:
        #    mut_type = self.mutations[-1].mut_type
        #else:
        #    mut_type = 'n'
        mut_type = 'n'

        if self.size > 0:
            self_node = [(self.pop_count(),
                          "pr-{}-{}{}".format(str(self.prolif_rate),
                                              mut_type, idnt))]
            #self.idnt = idnt + str(self.depth)
            #if len(self.nodes) == 0: #deepest point in tree
                #self_node = [

        if not self.nodes:
            return self_node
        else:
            for node in self.nodes:
                # TODO fix this - it's currently broken.
                #if node.mutations:
                #    mut_type = node.mutations[-1].mut_type
                #else:
                #    mut_type = 'n'
                mut_type = 'n'
                new_idnt = str(self.nodes.index(node)) + mut_type + idnt
                blank_ctree += node.freq_to_list(new_idnt)
            return blank_ctree + self_node

    def get_clone_attrs_as_list(self, attr_names, inc_dead_clones=False):
        """
        Get specified attributes of clone, and its children, as a flat list.
        """
        if isinstance(attr_names, basestring):
            attr_names = [attr_names]

        # get attributes of this particular clone
        attr_list = []
        if not self.is_dead() or inc_dead_clones:
            for name in attr_names:
                try:
                    attr_val = getattr(self, name)
                except AttributeError:
                    raise
                attr_list.append(attr_val)

        if len(attr_list) > 1:
            # multiple attributes have been recorded;
            # need to condense to a single list element
            attr_list = [tuple(attr_list)]

        # get attributes of children, if any
        if self.has_children():
            for node in self.nodes:
                # note the need to use 'splat' operator to unpack attr_names
                attr_list += node.get_clone_attrs_as_list(attr_names,
                                                          inc_dead_clones)
        return attr_list

    def is_dead(self):
        """Determine if this clone is dead, i.e. has no cells."""
        return self.size <= 0

    def has_children(self):
        """Determine if this clone has children in the subpopulation tree."""
        return bool(self.nodes)

    def is_dead_end(self):
        """Determine if this clone is a dead end."""
        return self.is_dead() and not self.has_children()

    def to_JSON(self):
        return json.dumps(self, default=subpop_to_JSON,
                          sort_keys=True, indent=4)



def subpop_to_JSON(obj):
    if isinstance(obj, Subpopulation):
        d = {}
        d['__class__'] = 'Subpopulation'
        d['__value__'] = obj.__dict__
        return d
    else:
        raise TypeError(repr(obj) + ' is not JSON serializable')


def safe_binomial_sample(num, prob):
    """Sample from the binomial distribution, validating params first."""
    if num < 0 or prob < 0:
        return 0
    else:
        return np.random.binomial(num, prob)
