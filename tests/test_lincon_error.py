"""
Tests that pyOptSparse raises an error when a user-defined obj/con function returns a linear constraint value
(which should not because linear constraint is defined exclusively by jac and bounds)
"""

# Standard Python modules
import unittest

# First party modules
from pyoptsparse import SLSQP, Optimization


def objfunc(xdict):
    """Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3"""
    x = xdict["x"]
    funcs = {}

    funcs["obj"] = x**2
    funcs["con"] = x - 1  # falsely return a linear constraint value

    fail = False
    return funcs, fail


class TestLinearConstraintCheck(unittest.TestCase):
    def test(self):
        # define an optimization problem with a linear constraint
        optProb = Optimization("test", objfunc)
        optProb.addVarGroup("x", 1, value=1)
        optProb.addObj("obj")
        optProb.addConGroup("con", 1, lower=1.0, linear=True, wrt=["x"], jac={"x": [1.0]})

        opt = SLSQP()
        with self.assertRaises(ValueError) as context:
            opt(optProb, sens="FD")

        # check if we get the expected error message
        err_msg = (
            "Value for linear constraint returned from user obj function. Linear constraints "
            + "are evaluated internally and should not be returned from the user's function."
        )
        self.assertEqual(err_msg, str(context.exception))


if __name__ == "__main__":
    unittest.main()
