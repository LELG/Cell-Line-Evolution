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
    def __init__(self, opt, mut_id=None):
        """Create new mutation"""
        if not mut_id:
            self.mut_id = id(self)
        self.prolif_rate_effect = get_prolif_rate_mutn(opt)

def get_prolif_rate_mutn(opt):
    """Generate a proliferation rate mutation effect."""
    mutn_magnitude = opt.sample_from_beta_dist() * opt.scale
    sample_for_mut_type = random.random() - opt.prob_del_mut
    mutn_valence = math.copysign(1, sample_for_mut_type)
    return mutn_magnitude * mutn_valence

