"""
Module for class representing mutations,
and associated functions.
"""
import math
import numpy as np
import random

class Mutation(object):
    """
    Class to represent mutations.
    """
    def __init__(self, opt, all_muts, mut_id=None):
        """Create new mutation"""
        if not mut_id:
            self.mut_id = id(self)
        # get proliferation rate effect and mutation type
        self.prolif_rate_effect = get_prolif_rate_mutn(opt)
        if self.prolif_rate_effect == 0.0:
            self.mut_type = 'n'
        elif self.prolif_rate_effect > 0.0:
            self.mut_type = 'b'
        else:
            self.mut_type = 'd'
        # get mutation rate effect - this does not
        # affect mutation type
        self.mut_rate_effect = get_mut_rate_mutn(opt)

        self.original_clone = None
        self.resist_strength = None

        # assume that all_muts is a dictionary with
        # keys for each mutation type
        try:
            all_muts[self.mut_type].append(self)
        except:
            raise

    def switch_mutn_type(self, all_muts, new_mut_type):
        """Change the mutation type of this mutation."""
        try:
            all_muts[self.mut_type].remove(self)
        except ValueError:
            raise ValueError("Mutn not in all_muts[{}]".format(self.mut_type))
        try:
            all_muts[new_mut_type].append(self)
        except KeyError:
            raise
        self.mut_type = new_mut_type

    def become_resistant(self, all_muts, resist_strength):
        """Make this mutation a resistance mutation."""
        # switch which sublist this mutation appears in
        self.original_clone.switch_mutn_type(self, 'r')
        self.switch_mutn_type(all_muts, 'r')
        self.resist_strength = resist_strength
        self.original_clone.become_resistant(self.resist_strength)


def get_prolif_rate_mutn(opt):
    """Generate a proliferation rate mutation effect."""
    scale_factor = opt.scale * opt.pro
    prolif_rate_effect = get_mutn_effect(opt.get_beta_dist_sample,
                                         scale_factor, opt.prob_mut_pos,
                                         opt.prob_mut_neg)
    return prolif_rate_effect


def get_mut_rate_mutn(opt):
    """Generate a mutation rate mutation effect."""
    scale_factor = opt.mscale * opt.mut
    mut_rate_effect = get_mutn_effect(opt.get_beta_dist_sample,
                                      scale_factor, opt.prob_inc_mut,
                                      opt.prob_dec_mut)
    return mut_rate_effect


def get_mutn_effect(get_effect_size, scale_factor, prob_pos, prob_neg):
    """Get a mutation effect size and type.

    Note that `get_effect_size` must be a function
    which samples from a suitable probability distribution.
    """
    # get a random effect size from a suitable probability
    # distribution, and scale as necessary
    mutn_magnitude = get_effect_size() * scale_factor

    # random.random generates a uniform pseudo-random float
    # between 0 and 1; by subtracting the probability
    # of a negative mutation from this float, we get
    # a float with `prob_neg` chance of being < 0
    # and `1 - prob_neg` chance of being > 0
    # (which could be beneficial or neutral)
    mut_type_from_sign = random.random() - prob_neg

    # allow for the possibility that there is a non-zero
    # chance of a 'strictly' neutral mutation
    prob_neutral_mut = (1 - prob_pos - prob_neg)
    if 0.0 <= mut_type_from_sign < prob_neutral_mut:
        mutn_magnitude = 0.0

    # now we copy the sign of this random float (our
    # mutation 'type') to the effect size we got from
    # the random distribution
    mutn_effect = math.copysign(mutn_magnitude, mut_type_from_sign)

    return mutn_effect


def generate_resistance(all_mutations, tumoursize, deterministic_num_r_muts, resist_strength=1.0):
    """Trigger the creation of resistance mutations."""
    # create flat list of deleterious/neutral mutations
    del_neutr_mutns = all_mutations['d'] + all_mutations['n']

    if deterministic_num_r_muts >= 0:
        num_resist_mutns = deterministic_num_r_muts
    else:
        num_resist_mutns = get_rand_num_r_mutns(len(del_neutr_mutns), tumoursize)

    # random.sample() returns a list of num_resist_mutns mutations,
    # randomly selected from del_neutr_mutns (with uniform likelihood)
    resistance_mutns = random.sample(del_neutr_mutns, num_resist_mutns)

    for mutn in resistance_mutns:
        mutn.become_resistant(all_mutations, resist_strength)

    print("{} resistance mutations generated.".format(num_resist_mutns))


def get_rand_num_r_mutns(total_mutns, tumoursize, min_resistant_pop_size=1e6):
    """Get a random number of resistance mutations for a given population."""
    prob_single_resist_mutn = tumoursize / float(min_resistant_pop_size)
    prob_resistance_mutn = prob_single_resist_mutn / float(total_mutns)
    num_resist_mutns = np.random.binomial(total_mutns, prob_resistance_mutn)
    return num_resist_mutns
