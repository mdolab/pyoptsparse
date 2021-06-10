# External modules
import numpy as np


class Objective(object):
    def __init__(self, name, scale=1.0):
        """
        This class holds the representation of a pyOptSparse objective.

        Parameters
        ----------
        name : str
            Name of this objective

        scale : float
            Scaling factor for objective. This does not change the actual
            optimization problem, but may be used to give a more
            human-meaningful value

        See Also
        --------
        Optimization.addObj : for the full documentation
        """
        self.name = name
        self.value = 0.0
        self.scale = scale

    def __str__(self):
        """
        Structured Print of Objective
        """
        res = "        Name        Value\n"
        res += "	 " + str(self.name).center(9)
        res += "%12g\n" % (np.real(self.value))

        return res
