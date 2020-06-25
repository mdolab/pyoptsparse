#!/usr/bin/env python
import copy
from .pyOpt_optimization import Optimization


class Solution(Optimization):
    def __init__(self, optProb, optTime, optInform):
        """
        This class is used to describe the solution of an optimization
        problem. This class inherits from Optimization which enables a
        solution to be used as an input to a subsequent optimization problem.

        Parameters
        ----------
        optProb : Optimization problem class
            Optimization problem used to create solution

        optTime : float
            Time required for the optimization

        optInform : int
            The inform code from the optimization.
        """

        Optimization.__init__(self, optProb.name, None)

        # Copy over the variables, constraints, and objectives
        self.variables = copy.deepcopy(optProb.variables)
        self.constraints = copy.deepcopy(optProb.constraints)
        self.objectives = copy.deepcopy(optProb.objectives)
        self.optTime = optTime
        self.optInform = optInform

    def __str__(self):
        """
        Print Structured Solution
        """
        text0 = Optimization.__str__(self)
        text1 = ""
        lines = text0.split("\n")
        lines[1] = lines[1][len("Optimization Problem -- ") :]
        for i in range(5):
            text1 += lines[i] + "\n"

        text1 += "\n    Solution: \n"
        text1 += ("-" * 80) + "\n"
        text1 += "    Total Time: %25.4f\n" % self.optTime
        text1 += "       User Objective Time :   %10.4f\n" % self.userObjTime
        text1 += "       User Sensitivity Time : %10.4f\n" % self.userSensTime
        if hasattr(self, "userJProdTime"):
            text1 += "       User J-Product Time :   %10.4f\n" % self.userJProdTime
            text1 += "       User J^T-Product Time : %10.4f\n" % self.userJTProdTime
        text1 += "       Interface Time :        %10.4f\n" % self.interfaceTime
        text1 += "       Opt Solver Time:        %10.4f\n" % self.optCodeTime
        text1 += "    Calls to Objective Function : %7d\n" % self.userObjCalls
        text1 += "    Calls to Sens Function :      %7d\n" % self.userSensCalls
        if hasattr(self, "userJProdCalls"):
            text1 += "    Calls to JProd Function :     %7d\n" % self.userJProdCalls
            text1 += "    Calls to JTProd Function :    %7d\n" % self.userJTProdCalls

        for i in range(5, len(lines)):
            text1 += lines[i] + "\n"

        text1 += ("-" * 80) + "\n"

        return text1
