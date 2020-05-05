#!/usr/bin/env python
"""
pyOpt_objective

Holds the representation of a pyOptSparse objective
"""
import numpy as np

# =============================================================================
# Objective Class
# =============================================================================
class Objective(object):

    """
    Optimization Objective Class

    Parameters
    ----------
    name : str
        Name of this objective

    scale : float
        Scaling factor for objective. This does not change the actual
        optimization problem, but may be used to give a more
        human-meaningful value
        """

    def __init__(self, name, scale=1.0):
        self.name = name
        self.value = 0.0
        self.optimum = 0.0
        self.scale = scale

    def __str__(self):
        """
        Structured Print of Objective
        """
        res = "        Name        Value        Optimum\n"
        res += "	 " + str(self.name).center(9)
        res += "%12g  %12g\n" % (np.real(self.value), np.real(self.optimum))

        return res
