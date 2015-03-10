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
    def __init__(self, opt, subpop, all_muts, mut_id=None):
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

        self.original_clone = subpop
        # assume that all_muts is a dictionary with
        # keys for each mutation type
        try:
            all_muts[self.mut_type].append(self)
        except:
            raise

    def become_resistant(self):
        # get a strong beneficial mutation effect
        # TODO implement me
        # signal to clone that something has changed
        self.original_clone.recalculate_fitness()


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


def generate_resistance(all_mutations, tumoursize):
    """Trigger the creation of resistance mutations."""
    # TODO This will not work now that all_mutations is a dict
    total_mutns = len(all_mutations)
    num_resist_mutns = get_num_resist_mutations(total_mutns, tumoursize)
    # TODO need to change this so that the list is a flattened
    # TODO list of neutral/del mutations from all_mutations dict.
    # random.sample selects num_resist_mutns from the list
    # of all mutations (with uniform likelihood)
    resistance_mutns = random.sample(all_mutations, num_resist_mutns)
    for mutn in resistance_mutns:
        mutn.become_resistant()
    return resistance_mutns


def get_num_resist_mutations(total_mutns, tumoursize, min_resistant_pop_size=1e6):
    """Get a random number of resistance mutations for a given population."""
    prob_single_resist_mutn = tumoursize / float(min_resistant_pop_size)
    prob_resistance_mutn = prob_single_resist_mutn / float(total_mutns)
    num_resist_mutns = np.random.binomial(total_mutns, prob_resistance_mutn)
    return num_resist_mutns
