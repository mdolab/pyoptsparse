# Local modules
from .pyOpt_utils import INFINITY


class Variable:
    def __init__(
        self,
        name: str,
        varType: str,
        value,
        lower,
        upper,
        scale,
        offset,
        scalar=False,
        choices=[],
    ):
        """
        This class holds the representation of a single pyOptSparse variable

        See Also
        --------
        pyoptsparse.pyOpt_optimization.Optimization.addVarGroup : for the full documentation
        """
        self.name = name
        self.type = varType
        self.scalar = scalar
        self.choices = choices
        if self.type == "c":
            if lower is None:
                self.lower = -INFINITY
            else:
                self.lower = (lower - offset) * scale

            if upper is None:
                self.upper = INFINITY
            else:
                self.upper = (upper - offset) * scale

            self.value = (value - offset) * scale
            self.scale = scale
            self.offset = offset
        elif self.type == "i":
            self.value = int(value)
            self.lower = lower
            self.upper = upper
            self.scale = scale
            self.offset = offset
        elif self.type == "d":
            if len(choices) == 0:
                raise ValueError("A discrete variable requires to input an array of choices.")
            self.value = self.choices[int(value)]
            self.lower = 0
            self.upper = len(self.choices)
            self.scale = scale

    def __eq__(self, other):
        """
        Compare two variable objects
        """
        if (
            self.name == other.name
            and self.type == other.type
            and self.scalar == other.scalar
            and self.upper == other.upper
            and self.lower == other.lower
            and self.choices == other.choices
        ):
            return True
        else:
            return False

    def __str__(self) -> str:
        """
        Print Structured List of Variable
        """

        res = "Name     Type       Value       "
        res += "Lower Bound  Upper Bound\n"

        if self.type == "d":
            res += "	 "
            res += str(self.name).center(15)
            res += f"{self.type:>25}{self.choices[int(self.value)]:20f} {min(self.choices):14.2e} {max(self.choices):12.2e}\n"
        else:
            lower = self.lower
            upper = self.upper
            if self.lower is None:
                lower = -INFINITY
            if self.upper is None:
                upper = INFINITY

            res += "	 "
            res += str(self.name).center(9)
            res += f"{self.type:>5}	{self.value:14f} {lower:14.2e} {upper:12.2e} \n"

        return res
