"""
pyEgor - A pyOptSparse interface to Egor from egobox.
"""

# Standard Python modules
import datetime
import inspect
import time

# External modules
import numpy as np

# Local modules
from ..pyOpt_solution import SolutionInform
from ..pyOpt_optimizer import Optimizer
from ..pyOpt_utils import import_module

# import the Python module
egobox = import_module("egobox")


class Egor(Optimizer):
    """
    Egor Optimizer Class - wrapper for the Egor global optimization algorithm from egobox.
    """

    def __init__(self, raiseError=True, options=None):
        if options is None:
            options = {}
        name = "Egor"
        category = "Global Optimizer"
        defOpts = self._getDefaultOptions()
        informs = self._getInforms()
        super().__init__(name, category, defaultOptions=defOpts, informs=informs, options=options)

        if isinstance(egobox, Exception) and raiseError:
            raise egobox

    @staticmethod
    def _getInforms():
        informs = {
            1: "Reached maximum number of iterations",
            2: "Reached target cost function value",
            3: "Algorithm manually interrupted with SIGINT (Ctrl+C), SIGTERM or SIGHUP",
            4: "Algorithm peek at the same point twice. We consider it is converged.",
            5: "Timeout reached",
            6: "Solver unexpected exit. See logs for details.",
        }
        return informs

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            "gp_config": [dict, dict()],  # GpConfig as a dict used by Egor for surrogate model configuration
            "cstr_tol": [list, []],
            "n_start": [int, 20],
            "n_doe": [int, 0],
            "doe": [list, [[]]],
            "infill_strategy": [int, 4],  # default to LOG_EI
            "cstr_infill": [bool, False],
            "cstr_strategy": [int, 1],  # default to MC
            "qei_config": [dict, dict()],
            "infill_optimizer": [int, 1],  # default to COBYLA
            "trego": [dict, dict()],
            "coego_n_coop": [int, 0],
            "target": [float, -1e12],
            "outdir": [str, ""],
            "warm_start": [bool, False],
            "hot_start": [bool, False],
            "failsafe_strategy": [int, 1],  # default to REJECTION
            "seed": [int, -1],
            "verbose": [int, 0],  # level of verbosity, 0 = error, 1 = warn, 2 = info, 3 = debug
            "max_iters": [int, 20],
            "run_info": [dict, dict()],
            "timeout": [float, -1.0],
            "fcstrs": [list, []],
            "fcstr_specs": [list, []],
        }
        return defOpts

    def __call__(self, optProb, storeHistory=None, hotStart=None, **kwargs):
        """
        Solve optimization problem using Egor from egobox.

        Notes
        -----
        The kwargs are present for compatibility with other optimizers.
        Any sensitivity settings are ignored because Egor is derivative-free.
        """
        self.startTime = time.time()

        # Save the optimization problem and finalize constraint Jacobian
        self.optProb = optProb
        self.optProb.finalize()

        # Egor currently supports a single objective in this wrapper.
        if len(self.optProb.objectives) != 1:
            raise ValueError("Egor wrapper currently supports single-objective problems only.")

        # Set history/hotstart/coldstart
        self._setHistory(storeHistory, hotStart)
        self._setInitialCacheValues()

        if len(optProb.constraints) == 0:
            self.unconstrained = True

        blx, bux, xs = self._assembleContinuousVariables()
        xs = np.maximum(xs, blx)
        xs = np.minimum(xs, bux)
        n = len(xs)

        if np.any(~np.isfinite(blx)) or np.any(~np.isfinite(bux)):
            raise ValueError("Egor requires finite lower and upper bounds for all design variables.")

        if self.unconstrained:
            m = 0
        else:
            indices, blc, buc, fact = self.optProb.getOrdering(["ne", "le", "ni", "li"], oneSided=True, noEquality=True)
            m = len(indices)
            self.optProb.jacIndices = indices
            self.optProb.fact = fact
            self.optProb.offset = buc

        if self.optProb.comm.rank == 0:
            opt = self.getOption

            # Build x specifications from pyOptSparse bounds.
            xspecs = [egobox.XSpec(egobox.XType.FLOAT, [float(blx[i]), float(bux[i])]) for i in range(n)]

            gp_config = opt("gp_config")
            if gp_config is None:
                gp_config = egobox.GpConfig()

            infill_strategy = opt("infill_strategy")
            if infill_strategy is None:
                infill_strategy = egobox.InfillStrategy.LOG_EI

            cstr_strategy = opt("cstr_strategy")
            if cstr_strategy is None:
                cstr_strategy = egobox.ConstraintStrategy.MC

            qei_config = opt("qei_config")
            if qei_config is None:
                qei_config = egobox.QEiConfig()

            infill_optimizer = opt("infill_optimizer")
            if infill_optimizer is None:
                infill_optimizer = egobox.InfillOptimizer.COBYLA

            failsafe_strategy = opt("failsafe_strategy")
            if failsafe_strategy is None:
                failsafe_strategy = egobox.FailsafeStrategy.REJECTION

            fcstrs_opt = opt("fcstrs")
            n_cstr = m
            fcstr_specs = opt("fcstr_specs")

            n_fcstrs = 0 if fcstrs_opt is None else len(fcstrs_opt)
            if fcstr_specs is not None and len(fcstr_specs) not in (0, n_fcstrs):
                raise ValueError(
                    "Option 'fcstr_specs' length must be zero or match the number of function constraints."
                )

            ctor_supported = set(inspect.signature(egobox.Egor).parameters.keys())
            supports_ctor_verbose = "verbose" in ctor_supported

            if opt("verbose") is not None and not supports_ctor_verbose:
                raise ValueError("Installed egobox version does not support constructor option 'verbose'.")

            ctor_kwargs = {
                "gp_config": gp_config,
                "n_cstr": n_cstr,
                "cstr_tol": opt("cstr_tol") if len(opt("cstr_tol")) > 0 else None,
                "n_start": opt("n_start"),
                "n_doe": opt("n_doe"),
                "doe": np.array(opt("doe")) if np.array(opt("doe")).size > 0 else None,
                "infill_strategy": infill_strategy,
                "cstr_infill": opt("cstr_infill"),
                "cstr_strategy": cstr_strategy,
                "qei_config": qei_config,
                "infill_optimizer": infill_optimizer,
                "trego": opt("trego") if opt("trego") else None,
                "coego_n_coop": opt("coego_n_coop"),
                "target": float(opt("target")),
                "failsafe_strategy": failsafe_strategy,
            }
            ctor_kwargs = {k: v for k, v in ctor_kwargs.items() if k in ctor_supported}
            solver = egobox.Egor(xspecs, **ctor_kwargs)

            def fun(x):
                x_eval = np.atleast_2d(np.asarray(x, dtype=float))
                ncols = 1 + m
                y = np.zeros((x_eval.shape[0], ncols), dtype=float)
                for i in range(x_eval.shape[0]):
                    xi = np.clip(x_eval[i], blx, bux)
                    fobj, fcon, fail = self._masterFunc(xi, ["fobj", "fcon"])
                    if fail:
                        y[i, :] = np.nan
                        continue
                    y[i, 0] = float(np.atleast_1d(fobj)[0])
                    if m > 0:
                        y[i, 1:] = np.asarray(fcon, dtype=float)
                return y

            fcstrs = [] if fcstrs_opt is None else list(fcstrs_opt)

            minimize_kwargs = {
                "fcstrs": fcstrs,
                "fcstr_specs": [] if fcstr_specs is None else fcstr_specs,
                "max_iters": opt("max_iters"),
                "run_info": opt("run_info"),
                "outdir": opt("outdir") if opt("outdir") != "" else None,
                "warm_start": opt("warm_start"),
                "hot_start": True if opt("hot_start") else False,
                "seed": opt("seed") if opt("seed") >= 0 else None,
                "timeout": float(opt("timeout")) if opt("timeout") > 0 else None,
                "verbose": opt("verbose"),
            }

            t0 = time.time()
            egor_result = solver.minimize(fun, **minimize_kwargs)
            optTime = time.time() - t0

            if self.storeHistory:
                self.metadata["endTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.metadata["optTime"] = optTime
                self.hist.writeData("metadata", self.metadata)
                self.hist.close()

            # Broadcast a -1 to indicate optimization has finished
            self.optProb.comm.bcast(-1, root=0)

            # Optimizer has no standardized exit code mapping in this wrapper.
            sol_inform = SolutionInform.from_informs(self.informs, int(egor_result.status.exit))

            result = egor_result.result

            xstar = np.asarray(result.x_opt, dtype=float).reshape(-1)
            ystar = np.asarray(result.y_opt, dtype=float).reshape(-1)
            fstar = float(ystar[0])

            # Create the optimization solution
            sol = self._createSolution(optTime, sol_inform, fstar, xstar)
        else:
            self._waitLoop()
            sol = None

        sol = self._communicateSolution(sol)
        return sol
