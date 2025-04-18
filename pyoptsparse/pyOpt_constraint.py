# Standard Python modules
from collections import OrderedDict
import copy
from typing import Dict, Iterable, List, Optional, Union

# External modules
import numpy as np

# Local modules
from .pyOpt_error import pyOptSparseWarning
from .pyOpt_types import Dict1DType
from .pyOpt_utils import INFINITY, _broadcast_to_array, convertToCOO


class Constraint:
    def __init__(
        self,
        name: str,
        nCon: int,
        linear: bool,
        wrt: Union[None, str, Iterable[str]],
        jac: Dict1DType,
        lower,
        upper,
        scale,
    ):
        """
        This class holds the representation of a pyOptSparse constraint group

        See Also
        --------
        pyoptsparse.pyOpt_optimization.Optimization.addConGroup : for the full documentation
        """
        self.name = name
        self.ncon = nCon
        self.linear = linear
        self.wrt = wrt
        self.jac = jac
        self.partialReturnOk: Optional[bool] = None
        self.scale = scale
        self.rs: Optional[int] = None
        self.re: Optional[int] = None
        # Before we can do the processing below we need to have lower
        # and upper arguments expanded:

        lower = _broadcast_to_array("lower", lower, nCon, allow_none=True)
        upper = _broadcast_to_array("upper", upper, nCon, allow_none=True)
        scale = _broadcast_to_array("scale", scale, nCon)

        # Save lower and upper...they are only used for printing however
        self.lower = lower
        self.upper = upper
        # The current value of the constraint (for printing purposes)
        self.value = np.zeros(self.ncon)

        # Now we determine what kind of constraint this is:
        # 1. An equality constraint
        # 2. A upper bound on a 1-sided constraint
        # 3. A lower bound on a 1-sided constraint
        # 4. Lower and Upper bounds on 2-sided constraint
        # 5. No lower or upper bounds. Typically will only be used for
        # dummy constraint on an unconstrained problem.

        # The first 3, will give a single "constraint" in all
        # optimizers. Some optimizers can only do 1-sided constraints
        # so type 4 and 5 must be split into two separate constraints
        # automatically.

        # This keeps track of the equality constraints:
        equalityConstraints: Dict[str, List] = {"value": [], "ind": [], "fact": []}

        # All (inequality) constraints get added to
        # "twoSidedConstraints". This will be used in optimizers that
        # can do two-sided constraints properly
        twoSidedConstraints: Dict[str, List] = {"lower": [], "upper": [], "ind": [], "fact": []}

        # All (inequality) constraints are also added to
        # "oneSidedConstraints". These are processed such that the
        # lower bound is ALWAYS -INFINITY such that: con <= upper For
        # optimizers that need things <= zero, this can be processed
        # with a (-value) offset. One sided constraints need a fact
        # defined which is precisely 1.0 or -1.0. The -1.0 appears
        # when a greater-than-constraint is turned into a
        # less-than-constraint.
        oneSidedConstraints: Dict[str, List] = {"lower": [], "upper": [], "ind": [], "fact": []}

        for icon in range(self.ncon):
            # Check for equality constraint:
            if lower[icon] == upper[icon] and lower[icon] is not None:
                equalityConstraints["value"].append(lower[icon] * scale[icon])
                equalityConstraints["ind"].append(icon)
                equalityConstraints["fact"].append(1.0)

            # Two sided constraint:
            elif lower[icon] is not None and upper[icon] is not None:
                twoSidedConstraints["lower"].append(lower[icon] * scale[icon])
                twoSidedConstraints["upper"].append(upper[icon] * scale[icon])
                twoSidedConstraints["ind"].append(icon)
                twoSidedConstraints["fact"].append(1.0)

                # TWO sets of 1 sided constraints:
                oneSidedConstraints["lower"].append(-INFINITY)
                oneSidedConstraints["upper"].append(upper[icon] * scale[icon])
                oneSidedConstraints["ind"].append(icon)
                oneSidedConstraints["fact"].append(1.0)

                oneSidedConstraints["lower"].append(-INFINITY)
                oneSidedConstraints["upper"].append(-lower[icon] * scale[icon])
                oneSidedConstraints["ind"].append(icon)
                oneSidedConstraints["fact"].append(-1.0)

            # Upper bound only:
            elif upper[icon] is not None:
                twoSidedConstraints["lower"].append(-INFINITY)
                twoSidedConstraints["upper"].append(upper[icon] * scale[icon])
                twoSidedConstraints["ind"].append(icon)
                twoSidedConstraints["fact"].append(1.0)

                # Just one, 1-sided constraint
                oneSidedConstraints["lower"].append(-INFINITY)
                oneSidedConstraints["upper"].append(upper[icon] * scale[icon])
                oneSidedConstraints["ind"].append(icon)
                oneSidedConstraints["fact"].append(1.0)

            # Lower bound only:
            elif lower[icon] is not None:
                twoSidedConstraints["lower"].append(lower[icon] * scale[icon])
                twoSidedConstraints["upper"].append(INFINITY)
                twoSidedConstraints["ind"].append(icon)
                twoSidedConstraints["fact"].append(1.0)

                # Just one, 1-sided constraint
                oneSidedConstraints["lower"].append(-INFINITY)
                oneSidedConstraints["upper"].append(-lower[icon] * scale[icon])
                oneSidedConstraints["ind"].append(icon)
                oneSidedConstraints["fact"].append(-1.0)

            # Fully unconstrained!
            elif lower[icon] is None and upper[icon] is None:
                twoSidedConstraints["lower"].append(-INFINITY)
                twoSidedConstraints["upper"].append(INFINITY)
                twoSidedConstraints["ind"].append(icon)
                twoSidedConstraints["fact"].append(1.0)

                # Since this is just a dummy constraint, we only need
                # a single one....it can just be less than INFINITY
                oneSidedConstraints["lower"].append(-INFINITY)
                oneSidedConstraints["upper"].append(INFINITY)
                oneSidedConstraints["ind"].append(icon)
                oneSidedConstraints["fact"].append(1.0)
            # end if (con type)
        # end for (con loop)

        # Convert the stuff to arrays:
        oneSidedConstraints["ind"] = np.array(oneSidedConstraints["ind"], "intc")
        twoSidedConstraints["ind"] = np.array(twoSidedConstraints["ind"], "intc")
        equalityConstraints["ind"] = np.array(equalityConstraints["ind"], "intc")

        oneSidedConstraints["fact"] = np.array(oneSidedConstraints["fact"])
        twoSidedConstraints["fact"] = np.array(twoSidedConstraints["fact"])
        equalityConstraints["fact"] = np.array(equalityConstraints["fact"])

        equalityConstraints["value"] = np.array(equalityConstraints["value"])

        # Now save this information:
        self.equalityConstraints = equalityConstraints
        self.oneSidedConstraints = oneSidedConstraints
        self.twoSidedConstraints = twoSidedConstraints

    def finalize(self, variables: OrderedDict, dvOffset, index: int):
        """
        After the design variables have been finalized and the order
        is known we can check the constraint for consistency.

        Parameters
        ----------
        variables : OrderedDict
            The pyOpt variable list after they have been finalized.

        dvOffset : dict
            Design variable offsets from pyOpt_optimization

        index : int
            The starting index of this constraint in natural order

        Warnings
        --------
            This function should not need to be called by the user
        """

        # Set the row start and end
        self.rs = index
        self.re = index + self.ncon

        # First check if 'wrt' is supplied...if not we just take all
        # the dvGroups
        if self.wrt is None:
            self.wrt = list(variables.keys())
        else:
            # Sanitize the wrt input:
            if isinstance(self.wrt, str):
                self.wrt = [self.wrt.lower()]
            else:
                try:
                    self.wrt = list(self.wrt)
                except Exception:
                    raise TypeError(f"The 'wrt' argument to constraint '{self.name}' must be an iterable list")

            # We allow 'None' to be in the list...they are null so
            # just pop them out:
            self.wrt = [dvGroup for dvGroup in self.wrt if dvGroup is not None]

            # Now, make sure that each dvGroup the user supplied list
            # *actually* are variables
            for dvGroup in self.wrt:
                if dvGroup not in variables:
                    raise KeyError(
                        f"The supplied dvGroup '{dvGroup}' in 'wrt' for the {self.name} constraint, does not exist. "
                        + "It must be added with a call to addVar() or addVarGroup()."
                    )

            # Check for duplicates in wrt
            wrt_uniq = list(set(self.wrt))
            if len(wrt_uniq) < len(self.wrt):
                duplicate_vars = list({x for x in self.wrt if self.wrt.count(x) > 1})
                pyOptSparseWarning(
                    f"The constraint {self.name} was created with duplicate variables in 'wrt'. "
                    + "The following duplicates were automatically removed: "
                )
                for var in duplicate_vars:
                    print(f"\t\t{var}")
                self.wrt = wrt_uniq

        # Last thing for wrt is to reorder them such that dvGroups are
        # in order. This way when the Jacobian is assembled in
        # processDerivatives() the coorindate matrix will in the right
        # order.
        dvStart = []
        for dvGroup in self.wrt:
            dvStart.append(dvOffset[dvGroup][0])

        # This sort wrt using the keys in dvOffset
        self.wrt = [x for (y, x) in sorted(zip(dvStart, self.wrt))]

        # Now we know which dvGroups this constraint will have a
        # derivative with respect to (i.e. what is in the wrt list)

        # Now, it is possible that Jacobians were given for none, some
        # or all the dvGroups defined in wrt.
        if self.jac is None:
            # If the constraint is linear we have to *Force* the user to
            # supply a constraint Jacobian for *each* of the values in
            # wrt. Otherwise, a matrix of zeros isn't meaningful for the
            # sparse constraints.

            if self.linear:
                raise ValueError(
                    "The 'jac' keyword to argument to addConGroup() must be supplied for a linear constraint. "
                    + f"The constraint in error is {self.name}."
                )

            # without any additional information about the Jacobian
            # structure, we must assume they are all dense.
            self.jac = {}
            for dvGroup in self.wrt:
                ss = dvOffset[dvGroup]
                ndvs = ss[1] - ss[0]
                self.jac[dvGroup] = convertToCOO(np.zeros((self.ncon, ndvs)))

            # Set a flag for the constraint object, that not returning
            # them all is ok.
            self.partialReturnOk = True

        else:
            # First sanitize input:
            if not isinstance(self.jac, dict):
                raise TypeError(
                    "The 'jac' keyword argument to addConGroup() must be a dictionary. "
                    + f"The constraint in error is {self.name}."
                )

            # Now loop over the set we *know* we need and see if any
            # are in jac. We will actually pop them out, and that way
            # if there is anything left at the end, we can tell the
            # user supplied information was unused.
            tmp = copy.deepcopy(self.jac)
            self.jac = {}
            for dvGroup in self.wrt:
                ss = dvOffset[dvGroup]
                ndvs = ss[1] - ss[0]

                try:
                    self.jac[dvGroup] = tmp.pop(dvGroup)
                except KeyError:
                    # No big deal, just make a dense component...and
                    # set to zero
                    self.jac[dvGroup] = convertToCOO(np.zeros((self.ncon, ndvs)))

                # Convert Now check that the supplied Jacobian to COO:
                self.jac[dvGroup] = convertToCOO(self.jac[dvGroup])

                # Generically check the shape:
                if self.jac[dvGroup]["shape"][0] != self.ncon or self.jac[dvGroup]["shape"][1] != ndvs:
                    raise ValueError(
                        f"The supplied Jacobian for dvGroup {dvGroup}' in constraint {self.name}, was the incorrect size. "
                        + f"Expecting a Jacobian of size ({self.ncon}, {ndvs}) but received a Jacobian of size "
                        + f"({self.jac[dvGroup]['shape'][0]}, {self.jac[dvGroup]['shape'][1]})."
                    )
            # end for (dvGroup)

            # If there is anything left in jac print a warning:
            for dvGroup in tmp:
                pyOptSparseWarning(
                    f"A Jacobian with dvGroup key of '{dvGroup}' was unused in constraint {self.name}. "
                    + "This will be ignored."
                )

            # Since this function *may* be called multiple times, only
            # set paritalReturnOk if it was the first pass:
            if self.partialReturnOk is None:
                # Finally partial returns NOT ok, since the user has
                # supplied information about the sparsity:
                self.partialReturnOk = False

        # end if (if Jac)

    def __str__(self) -> str:
        """
        Print Constraint
        """
        res = ""
        for i in range(self.ncon):
            lower = self.lower[i]
            upper = self.upper[i]
            value = self.value[i]
            if lower is None:
                lower = -INFINITY
            if upper is None:
                upper = INFINITY

            res += (
                "	 "
                + str(self.name).center(9)
                + f"	  i {np.real(lower):15.2e} <= {np.real(value):8f} <= {np.real(upper):8.2e}\n"
            )

        return res
