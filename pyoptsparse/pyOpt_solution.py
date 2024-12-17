# Standard Python modules
import copy

# External modules
import numpy as np

# Local modules
from .pyOpt_optimization import Optimization


class Solution(Optimization):
    def __init__(self, optProb, xStar, fStar, lambdaStar, optInform, info):
        """
        This class is used to describe the solution of an optimization
        problem. This class inherits from Optimization which enables a
        solution to be used as an input to a subsequent optimization problem.

        Parameters
        ----------
        optProb : Optimization problem class
            Optimization problem used to create solution

        xStar : dict
            The final design variables

        fStar : dict
            The final objective(s)

        lambdaStar : dict
            The final Lagrange multipliers

        optInform : int
            The inform code returned by the optimizer

        info : dict
            A dictionary containing timing and call counter info to be stored
            in the Solution object.
        """

        super().__init__(optProb.name, None)

        # Copy over the variables, constraints, and objectives
        self.variables = copy.deepcopy(optProb.variables)
        self.constraints = copy.deepcopy(optProb.constraints)
        self.objectives = copy.deepcopy(optProb.objectives)
        xopt = optProb._mapXtoOpt(optProb.processXtoVec(xStar))
        # Now set the x-values:
        i = 0
        for dvGroup in self.variables.keys():
            for var in self.variables[dvGroup]:
                var.value = xopt[i]
                i += 1

        # Now set the f-values
        if isinstance(fStar, np.ndarray) and len(fStar) == 1:
            self.objectives[list(self.objectives.keys())[0]].value = fStar.item()
            fStar = fStar.item()
        else:
            for f_name, f in self.objectives.items():
                f.value = fStar[f_name]

        self.optTime = info["optTime"]
        self.userObjTime = info["userObjTime"]
        self.userSensTime = info["userSensTime"]
        self.userObjCalls = info["userObjCalls"]
        self.userSensCalls = info["userSensCalls"]
        self.interfaceTime = info["interfaceTime"]
        self.optCodeTime = info["optCodeTime"]
        self.optInform = optInform
        self.fStar = fStar
        self.xStar = xStar
        self.lambdaStar = lambdaStar

    def __str__(self) -> str:
        """
        Print Structured Solution
        """
        text0 = self.summary_str(minimal_print=False, print_multipliers=True)
        text1 = ""
        lines = text0.split("\n")
        lines[1] = lines[1][len("Optimization Problem -- ") :]
        for i in range(5):
            text1 += lines[i] + "\n"

        text1 += "\n    Solution: \n"
        text1 += ("-" * 80) + "\n"
        text1 += f"    Total Time: {self.optTime:25.4f}\n"
        text1 += f"       User Objective Time :   {self.userObjTime:10.4f}\n"
        text1 += f"       User Sensitivity Time : {self.userSensTime:10.4f}\n"
        text1 += f"       Interface Time :        {self.interfaceTime:10.4f}\n"
        text1 += f"       Opt Solver Time:        {self.optCodeTime:10.4f}\n"
        text1 += f"    Calls to Objective Function : {self.userObjCalls:7}\n"
        text1 += f"    Calls to Sens Function :      {self.userSensCalls:7}\n"

        for i in range(5, len(lines)):
            text1 += lines[i] + "\n"

        inform_val = self.optInform["value"]
        inform_text = self.optInform["text"]
        text1 += "\n"
        text1 += "   Exit Status\n"
        text1 += "      Inform  Description\n"
        text1 += f"      {inform_val:>6}  {inform_text:<0}\n"

        text1 += ("-" * 80) + "\n"

        return text1
