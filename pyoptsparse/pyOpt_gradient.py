#!/usr/bin/env python
"""
pyOpt_gradient - A class that produce gradients using finite
difference or complex step.

Copyright (c) 2013-2014 by Dr. Gaetan Kenway
All rights reserved.

Developers
----------
- Dr. Gaetan Kenway (GKK)

History
-------
    v. 0.1  - Initial Class Creation (GKK)
"""
# =============================================================================
# External Python modules
# =============================================================================
import numpy
from .pyOpt_MPI import MPI
# =============================================================================
# Gradient Class
# =============================================================================
class Gradient(object):
    """
    Gradient class for automatically computing gradients with finite
    difference or complex step.

    Parameters
    ----------
    optProb : Optimization instance
        This is the complete description of the optimization problem.

    sensType : str
        'FD' for forward difference, 'CD' for central difference,
        'FDR' for forward difference with relative step size,
        'CDR' for central difference with relative step size,
        and 'CS' for complex step

    sensStep : number
        Step size to use for differencing

    sensMode : str
        Flag to compute gradients in parallel.
        """

    def __init__(self, optProb, sensType, sensStep=None, sensMode='',
                 comm=None):
        self.optProb = optProb
        self.sensType = sensType
        if sensStep is None:
            if self.sensType in ['fd', 'fdr']:
                self.sensStep = 1e-6
            elif self.sensType in ['cd', 'cdr']:
                self.sensStep = 1e-4
            else:
                self.sensStep = 1e-40j
        else:
            self.sensStep = sensStep
        self.sensMode = sensMode
        self.comm = comm

        # Now we can compute which dvs each process will need to
        # compute:
        ndvs = self.optProb.ndvs
        if self.sensMode == 'pgc' and self.comm:
            self.mydvs = list(range(self.comm.rank, ndvs, self.comm.size))
        else:
            self.mydvs = list(range(ndvs))

    def _eval_func(self, x):
        """internal method to call function and extract obj, con"""

        xCall = self.optProb.processX(x)
        # Call objective
        [funcs, fail] = self.optProb.objFun(xCall)

        # Process constraint in case they are in dict form
        self.optProb.evaluateLinearConstraints(x, funcs)
        fobj = self.optProb.processObjective(funcs, scaled=False)

        if self.sensType == 'cs':
            fcon = self.optProb.processConstraints(funcs, scaled=False,
                                                   dtype='D', natural=True)
        else:
            fcon = self.optProb.processConstraints(funcs, scaled=False,
                                                   natural=True)

        return fobj, fcon, fail

    def __call__(self, x, funcs):
        """
        We need to make this object "look" the same as a user supplied
        function handle. That way, the optimizers need not care how
        the gradients are **actually** calculated.

        Parameters
        ----------
        x : array
            Optimization variables from optimizer

        funcs : dict
            Dictionary of all the function values

        Returns
        -------
        gobj : 1D array
            The derivative of the objective with respect to the design
            variables

        gcon : 2D array
            The derivative of the constraints with respect to the design
            variables

        fail : bool
            Flag for failure. It currently always returns False
            """

        # Since this is *very* dumb loop over all the design
        # variables, it is easier to just loop over the x values as an
        # array. Furthermore, since the user **should** have
        # reasonably well scaled variables, the fixed step size should
        # have more meaning.

        # Generate final array sizes for the objective and constraint
        # gradients
        ndvs = self.optProb.ndvs
        ncon = self.optProb.nCon
        gobj = numpy.zeros(ndvs, 'd')
        gcon = numpy.zeros((ncon, ndvs), 'd')

        if self.sensMode == 'pgc':
            funcsBase = self.comm.bcast(funcs)
        else:
            funcsBase = funcs

        # We DO NOT want the constraints scaled here....the constraint
        # scaling will be taken into account when the derivatives are
        # processed as per normal.
        xBase = self.optProb.deProcessX(x)
        self.optProb.evaluateLinearConstraints(xBase, funcsBase)
        fconBase = self.optProb.processConstraints(
            funcsBase, scaled=False, dtype='D', natural=True)
        fobjBase = self.optProb.processObjective(funcsBase, scaled=False)

        # Convert to complex if necessary:
        if self.sensType == 'cs':
            xBase = xBase.astype('D')

        masterFail = False

        for i in self.mydvs:
            xph = xBase.copy()
            if self.sensType in ['fdr', 'cdr']:
                sensStep = max(abs(self.sensStep*xph[i]), self.sensStep)
            else:
                sensStep = self.sensStep
            xph[i] += sensStep

            fobj_ph, fcon_ph, fail = self._eval_func(xph)
            if fail:
                masterFail = True

            # forward difference
            if self.sensType in ['fd', 'fdr']:
                gobj[i]    = (fobj_ph - fobjBase)/sensStep
                gcon[:, i] = (fcon_ph - fconBase)/sensStep

            # central difference
            elif self.sensType in ['cd', 'cdr']:
                xmh = xph  # xph isn't used anymore so  just point to same location
                xmh[i] -= 2*sensStep

                fobj_mh, fcon_mh, fail = self._eval_func(xmh)
                if fail:
                    masterFail = True

                gobj[i]    = (fobj_ph - fobj_mh)/(2*sensStep)
                gcon[:, i] = (fcon_ph - fcon_mh)/(2*sensStep)

            # complex step
            else:
                gobj[i]    = numpy.imag(fobj_ph)/numpy.imag(sensStep)
                gcon[:, i] = numpy.imag(fcon_ph)/numpy.imag(sensStep)

        if self.sensMode == 'pgc':
            # We just mpi_reduce to the root with sum. This uses the
            # efficent numpy versions
            self.comm.Reduce(gobj.copy(), gobj, op=MPI.SUM, root=0)
            self.comm.Reduce(gcon.copy(), gcon, op=MPI.SUM, root=0)

        # Logically reduce (over the comm) if the fail if *ANY*
        # gradient calc failed:
        if self.comm is not None:
            masterFail = self.comm.allreduce(masterFail, op=MPI.LOR)

        # Finally, we have to convert everything **back** to a
        # dictionary so the rest of the code works:
        funcs = {}
        for objKey in self.optProb.objectives:
            funcs[objKey] = {}
            for dvGroup in self.optProb.variables:
                ss = self.optProb.dvOffset[dvGroup]
                funcs[objKey][dvGroup] = gobj[ss[0]:ss[1]]

        for conKey in self.optProb.constraints:
            con = self.optProb.constraints[conKey]
            funcs[conKey] = {}
            for dvGroup in self.optProb.variables:
                ss = self.optProb.dvOffset[dvGroup]
                funcs[conKey][dvGroup] = gcon[con.rs:con.re,ss[0]:ss[1]]

        return funcs, masterFail

