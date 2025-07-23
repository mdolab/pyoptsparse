"""
==============================================================================
ParOpt Test
==============================================================================
@File    :   test_paropt.py
@Date    :   2025/07/23
@Author  :   Alasdair Christison Gray
@Description :
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
import scipy.sparse as sp
import parameterized

# ==============================================================================
# Extension modules
# ==============================================================================
from pyoptsparse import Optimization, OPT
from testing_utils import OptTest

# We want to test the following functionality in ParOpt:
# - Dense vs Sparse constraint treatment
# - Trust region, interiot point, MMA algorithms
# - Constrained and unconstrained problems

test_params = []
for sparse in [True, False]:
    for constrained in [True, False]:
        for algorithm in ["tr", "ip", "mma"]:
            if not (sparse and algorithm == "tr"):
                test_params.append((sparse, constrained, algorithm))


def custom_name_func(testcase_func, param_num, param):
    sparseString = "Sparse" if param["sparse"] else "Dense"
    constrainedString = "Constrained" if param["constrained"] else "Unconstrained"
    algorithmString = param["algorithm"].upper()
    return f"{sparseString}_{constrainedString}_{algorithmString}"


@parameterized.parameterized_class(
    ("sparse", "constrained", "algorithm"),
    test_params,
    custom_name_func,
)
class TestParOpt(OptTest):
    def setUp(self):
        super().setUp()
        self.name = "Rosenbrock" if self.constrained else "Unconstrained Rosenbrock"
        self.N = 50
        self.cons = {"normCon"} if self.constrained else set()
        self.objs = {"obj"}
        self.DVs = {"x"}
        self.fStar = 214.380322 if self.constrained else 0.0
        self.xStar = {"x": np.ones(self.N) * 1 / np.sqrt(2.0) if self.constrained else np.ones(self.N)}
        self.optName = "ParOpt"

    def objfunc(self, xdict):
        x = xdict["x"]

        funcs = {}

        # for i in range(len(x) - 1):
        #     funcs["obj"] += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2
        funcs["obj"] = np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2)

        if self.constrained:
            # Nonlinear constraints, constraints are:
            # x0^2 + x1^2 <= 1
            # x1^2 + x2^2 <= 1
            # ...
            funcs["normCon"] = x[:-1] ** 2 + x[1:] ** 2

        fail = False
        return funcs, fail

    def sens(self, xdict, funcs):
        x = xdict["x"]
        funcsSens = {}
        grads = np.zeros(len(x))

        for i in range(len(x) - 1):
            grads[i] += 2 * (200 * x[i] ** 3 - 200 * x[i] * x[i + 1] + x[i] - 1)
            grads[i + 1] += 200 * (x[i + 1] - x[i] ** 2)

        funcsSens["obj"] = {"x": grads}

        if self.constrained:
            # Nonlinear constraints, constraints are:
            # x0^2 + x1^2 <= 1
            # x1^2 + x2^2 <= 1
            # ...
            dgdx = np.zeros((self.N - 1, self.N))
            for i in range(self.N - 1):
                dgdx[i, i] = 2 * x[i]
                dgdx[i, i + 1] = 2 * x[i + 1]
            dgdx = sp.coo_array(dgdx, shape=(self.N - 1, self.N))
            funcsSens["normCon"] = {"x": dgdx}

        fail = False
        return funcsSens, fail

    def setup_optProb(self):
        rng = np.random.default_rng(1234)

        optProb = Optimization("Rosenbrock Problem", self.objfunc)
        optProb.addObj("obj")
        optProb.addVarGroup("x", self.N, lower=0.5, upper=4.0, value=rng.uniform(0.5, 4.0, self.N))

        if self.constrained:
            # Linear constraints: x0 <= x1 <= x2 <= ... <= xN
            # x0 - x1 <= 0
            # x1 - x2 <= 0
            # ...
            rowInds = []
            colInds = []
            values = []
            for i in range(self.N - 1):
                rowInds += [i, i]
                colInds += [i, i + 1]
                values += [1.0, -1.0]
            jac = {"coo": [np.array(rowInds), np.array(colInds), np.array(values)], "shape": [self.N - 1, self.N]}
            optProb.addConGroup(
                "lincon",
                self.N - 1,
                upper=0.0,
                lower=-1.0,
                linear=True,
                jac={"x": jac},
            )

            # Nonlinear constraints have same sparsity structure as linear constraints
            optProb.addConGroup("normCon", nCon=self.N - 1, upper=1.0, jac={"x": jac})

        self.optProb = optProb

    def setup_optimizer(self, optOptions=None):
        options = {
            "algorithm": self.algorithm,
        }
        if optOptions is not None:
            options.update(optOptions)
        return OPT("ParOpt", sparse=self.sparse, options=options)

    def test_opt(self):
        self.setup_optProb()
        sol = self.optimize()
        self.assert_solution_allclose(sol, 1e-5)


if __name__ == "__main__":
    import unittest

    unittest.main()
    # from pyoptsparse import OPT

    # optProb, sens = generateOptProb(100, True)  # Example with 10 variables and constraints
    # # opt = OPT("ParOpt", sparse=False, options={"algorithm": "ip"})
    # opt = OPT("SNOPT")
    # sol = opt(optProb, sens=sens)
    # print("Solution:", sol)
