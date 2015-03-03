"""
Module for class representing mutations,
and associated functions.
"""
import math
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
        all_muts.append(self)


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

    # random.random generates a pseudo-random float
    # between 0 and 1; by subtracting the probability
    # of a negative mutation from this float, we get
    # a float with `prob_neg` chance of being negative
    # and `1 - prob_neg` chance of being positive
    #(which could mean beneficial, or neutral)
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
