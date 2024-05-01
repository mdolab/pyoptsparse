"""
alpso - Python Version of the Augmented Lagrangian Particle Swarm Optimizer

alpso if a global optimizer which solves problems of the form:

            min F(x)

    subject to: Gi(x)  = 0, i = 1(1)ME
                Gj(x) <= 0, j = ME+1(1)M
                xLB <= x <= xUB

To Do:
    - Migrate Inner Loop Printing Option
    - Add Other Inertia and Velocity Updates to Inner Loop
    - Fix Neighbourhood best from Lagrangian value
"""

# Standard Python modules
from math import floor
import os
import random
import time

# External modules
import numpy as np

# Local modules
from ..pyOpt_error import pyOptSparseWarning

# Misc Definitions
inf = 10.0e20  # define a value for infinity
eps = 1.0  # define a value for machine precision
while (eps / 2.0 + 1.0) > 1.0:
    eps /= 2.0

eps *= 2.0


# eps = math.ldexp(1,-52)


# ==============================================================================
# alpso function
# ==============================================================================
# fmt: off
def alpso(dimensions, constraints, neqcons, xtype, x0, xmin, xmax, swarmsize, nhn,
          nhm, maxOutIter, maxInnIter, minInnIter, stopCriteria, stopIters, etol,
          itol, rtol, atol, dtol, prtOutIter, prtInnIter, r0, vinit, vmax, c1, c2, w1, w2,
          ns, nf, vcrazy, fileout, filename, logfile, hstfile, rseed, scale, nhs, objfunc):
# fmt: on # noqa: E115
    """
    Python Version of the Augmented Lagrangian Particle Swarm Optimizer

    Documentation last updated:  April. 29, 2008 - Ruben E. Perez
    """

    #
    if x0.size > 0:
        if isinstance(x0, list):
            x0 = np.array(x0)
        elif not isinstance(x0, np.ndarray):
            pyOptSparseWarning(
                "Initial x must be either list or numpy.array, all initial positions randomly generated"
            )

    #
    if hstfile is not None:
        h_start = True
    else:
        h_start = False

    if logfile is not None:
        sto_hst = True
    else:
        sto_hst = False

    # Set random number seed
    rand = random.Random()
    if rseed == {}:
        rseed = time.time()

    rand.seed(rseed)

    #
    if filename == "":
        filename = "ALPSO.out"

    ofname = ""
    sfname = ""
    fntmp = filename.split(".")
    if len(fntmp) == 1:
        ofname += fntmp[0] + "_print.out"
        sfname += fntmp[0] + "_summary.out"
    else:
        if "/" not in fntmp[-1] and "\\" not in fntmp[-1]:
            ofname += filename[: filename.rfind(".")] + "_print." + fntmp[-1]
            sfname += filename[: filename.rfind(".")] + "_summary." + fntmp[-1]
        else:
            ofname += filename + "_print.out"
            sfname += filename + "_summary.out"

    header = ""
    header += " " * 37 + "======================\n"
    header += " " * 39 + " ALPSO 1.1 (Bulk)\n"
    header += " " * 37 + "======================\n\n"
    header += "Parameters:\n"
    header += "-" * 97 + "\n"
    if maxInnIter != minInnIter:
        diI = 1
    else:
        diI = 0

    if x0.size > 0:
        if len(x0.shape) == 1:
            nxi = 1
        else:
            nxi = x0.shape[0]

    else:
        nxi = 0
# fmt: off
    header += 'Swarmsize           :%9d' % swarmsize + '    MaxOuterIters     :%9d' % maxOutIter + '    Seed:%26.8f\n' % rseed
    header += 'Cognitive Parameter :%9.3f' % c1 + '    MaxInnerIters     :%9d' % maxInnIter + '    Scaling            :%11d\n' % scale
    header += 'Social Parameter    :%9.3f' % c2 + '    MinInnerIters     :%9d' % minInnIter + '    Stopping Criteria  :%11d\n' % stopCriteria
    header += 'Initial Weight      :%9.3f' % w1 + '    DynInnerIters     :%9d' % diI + '    Number of Failures :%11d\n' % ns
    header += 'Final Weight        :%9.3f' % w2 + '    StoppingIters     :%9d' % stopIters + '    Number of Successes:%11d\n\n' % nf

    header += 'Absolute Tolerance  : %1.2e' % atol + '    Number Initial Pos:%9d' % nxi + '    Neighbourhood Model:%11s\n' % nhm
    header += 'Relative Tolerance  : %1.2e' % rtol + '    Initial Velocity  :%9d' % vinit + '    Neighbourhood Size :%11d\n' % nhn
    header += 'Inequality Tolerance: %1.2e' % itol + '    Maximum Velocity  :%9d' % vmax + '    Selfless           :%11d\n' % nhs
    header += 'Equality Tolerance  : %1.2e' % etol + '    Craziness Velocity: %1.2e' % vcrazy + '    Fileout            :%11d\n' % fileout
    header += 'Global Distance     : %1.2e' % dtol + '    Initial Penalty   :%9.2f' % r0 + '    File Name          :%11s\n' % filename
    header += '-' * 97 + '\n\n'
# fmt: on
    if (fileout == 1) or (fileout == 3):
        if os.path.isfile(ofname):
            os.remove(ofname)

        ofile = open(ofname, "w")
        ofile.write(header)

    if (fileout == 2) or (fileout == 3):
        if os.path.isfile(sfname):
            os.remove(sfname)

        sfile = open(sfname, "w")
        sfile.write(header)

    #
    dt = 1.0
    vlimit = vmax
    vmax = np.ones(dimensions, float) * vmax
    if scale == 1:
        space_centre = np.zeros(dimensions, float)
        space_halflen = np.zeros(dimensions, float)
        for j in range(dimensions):
            space_centre[j] = (xmin[j] + xmax[j]) / 2.0
            space_halflen[j] = (xmax[j] - xmin[j]) / 2.0

        xmin = -np.ones(dimensions, float)
        xmax = np.ones(dimensions, float)
    else:
        for j in range(dimensions):
            vmax[j] = ((xmax[j] - xmin[j]) / 2.0) * vlimit

    # Initialize the positions and velocities for entire population
    x_k = np.zeros((swarmsize, dimensions), float)
    v_k = np.zeros((swarmsize, dimensions), float)
    discrete_i = []
    for i in range(swarmsize):
        for j in range(dimensions):
            x_k[i, j] = xmin[j] + rand.random() * (xmax[j] - xmin[j])
            if xtype[j] == 1:
                discrete_i.append(j)

            v_k[i, j] = (xmin[j] + rand.random() * (xmax[j] - xmin[j])) / dt

    if x0.size > 0:
        if len(x0.shape) == 1:
            if scale == 1:
                x_k[0, :] = (x0[:] - space_centre) / space_halflen
            else:
                x_k[0, :] = x0[:]

        else:
            if x0.shape[0] > swarmsize:
                pyOptSparseWarning(
                    "%d initial positions specified for %d particles, last %d positions ignored"
                    % (x0.shape[0], swarmsize, x0.shape[0] - swarmsize)
                )
                x0 = x0[0:swarmsize, :]

            for i in range(x0.shape[0]):
                if scale == 1:
                    x_k[i, :] = (x0[i, :] - space_centre) / space_halflen
                else:
                    x_k[i, :] = x0[i, :]

    # Initialize Augmented Lagrange
    f = np.zeros(swarmsize, float)
    L = np.zeros(swarmsize, float)
    g = np.zeros([swarmsize, constraints], float)
    g_old = np.zeros([swarmsize, constraints], float)
    rp = np.ones(constraints, float) * r0
    lambda_val = np.zeros(constraints, float)
    lambda_old = np.zeros(constraints, float)
    tau = np.zeros([swarmsize, constraints], float)
    tau_new = np.zeros(constraints, float)
    tau_old = np.zeros(constraints, float)
    nfevals = 0

    if h_start:
        [vals, hist_end] = hstfile.read([], ident=["obj", "con"])
        f = vals["obj"][0]
        g = vals["con"][0].reshape(g.shape)
    else:
        # Evaluate Objective Function
        if scale == 1:
            xtmp = (x_k * space_halflen) + space_centre
        else:
            xtmp = x_k

        for m in discrete_i:
            xtmp[:, m] = floor(xtmp[:, m] + 0.5)

        f, g = objfunc(xtmp)
        nfevals += swarmsize

    for i in range(swarmsize):

        # Augmented Lagrangian Value
        L[i] = f[i]
        if constraints > 0:

            # Equality Constraints
            for ell in range(neqcons):
                tau[i, ell] = g[i, ell]

            # Inequality Constraints
            for ell in range(neqcons, constraints):
                if rp[ell] != 0:
                    if g[i, ell] > -lambda_val[ell] / (2 * rp[ell]):
                        tau[i, ell] = g[i, ell]
                    else:
                        tau[i, ell] = -lambda_val[ell] / (2 * rp[ell])

                else:
                    tau[i, ell] = g[i, ell]

            #
            for ell in range(constraints):
                L[i] += lambda_val[ell] * tau[i, ell] + rp[ell] * tau[i, ell] ** 2

    # Initialize Particles Best
    best_x = np.zeros((swarmsize, dimensions))
    best_L = np.zeros(swarmsize, float)
    best_f = np.zeros(swarmsize, float)
    best_g = np.zeros([swarmsize, constraints], float)
    for i in range(swarmsize):
        for j in range(dimensions):
            best_x[i, j] = x_k[i, j]

        best_L[i] = L[i]
        best_f[i] = f[i]
        for ell in range(constraints):
            best_g[i, ell] = g[i, ell]

    # Initialize Swarm Best
    swarm_i = L.argmin()
    swarm_i_old = 0
    swarm_x = np.zeros(dimensions, float)
    for j in range(dimensions):
        swarm_x[j] = x_k[swarm_i, j]

    swarm_L = L[swarm_i]
    swarm_L_old = L[0]
    swarm_f = f[swarm_i]
    swarm_f_old = f[0]
    swarm_g = np.zeros(constraints, float)
    swarm_g_old = np.zeros(constraints, float)
    for ell in range(constraints):
        swarm_g[ell] = g[swarm_i, ell]
        swarm_g_old[ell] = g[0, ell]

    # Initialize Neighbourhood
    if (nhm == "dlring") or (nhm == "slring") or (nhm == "wheel") or (nhm == "spatial") or (nhm == "sfrac"):

        nhps = []
        nhbest_L = np.ones(swarmsize) * inf
        nhbest_f = np.zeros(swarmsize)
        nhbest_x = np.zeros((swarmsize, dimensions))
        nhbest_i = np.zeros(swarmsize)

        if nhm == "dlring":
            for i in range(swarmsize):
                nhps.append([])
                if nhs == 0:
                    nhps[i].append(i)

                for nb in range(1, (nhn / 2) + 1):
                    if i + nb >= swarmsize:
                        nhps[i].append(-1 + nb)
                    else:
                        nhps[i].append(i + nb)

                    if i - nb < 0:
                        nhps[i].append(swarmsize + i - nb)
                    else:
                        nhps[i].append(i - nb)

        elif nhm == "slring":
            for i in range(swarmsize):
                nhps.append([])
                if nhs == 0:
                    nhps[i].append(i)

                for nb in range(1, (nhn / 2) + 1):
                    if i + nb >= swarmsize:
                        nhps[i].append(-1 + nb)
                    else:
                        nhps[i].append(i + nb)

                    if i - (nb * 2) < 0:
                        nhps[i].append(swarmsize + i - (nb * 2))
                    else:
                        nhps[i].append(i - (nb * 2))

        elif nhm == "wheel":
            nhps.append([])
            nhps[0].append(0)
            for i in range(1, swarmsize):
                nhps.append([])
                nhps[i].append(i)
                nhps[i].append(0)
                nhps[0].append(i)

        elif nhm == "spatial":
            pdist = np.ones((swarmsize, swarmsize)) * inf
            for i in range(swarmsize):
                for i2 in range(i + 1, swarmsize):
                    pdist[i, i2] = np.linalg.norm(x_k[i2, :] - x_k[i, :])

                for i2 in range(i):
                    pdist[i, i2] = pdist[i2, i]

            for i in range(swarmsize):
                nhps.append([])
                for nb in range(nhn):
                    nhps[i].append(pdist[i, :].argmin())
                    pdist[i, nhps[i][nb]] = inf

                if nhs == 0:
                    nhps[i].append(i)

        elif nhm == "sfrac":
            pdist = np.zeros((swarmsize, swarmsize))
            d_max = np.zeros(swarmsize)
            frac = 0.6
            for i in range(swarmsize):
                for i2 in range(i + 1, swarmsize):
                    pdist[i, i2] = np.linalg.norm(x_k[i2, :] - x_k[i, :])

                for i2 in range(i):
                    pdist[i, i2] = pdist[i2, i]

            for i in range(swarmsize):
                nhps.append([])
                d_max[i] = pdist[i, :].max()
                for i2 in range(swarmsize):
                    if i == i2:
                        if nhs == 1:
                            pass
                        else:
                            nhps[i].append(i)

                    else:
                        if pdist[i, i2] / d_max[i] < frac:
                            nhps[i].append(i2)

        # Inizialize Neighbourhood Best
        for i in range(swarmsize):
            for nbp in nhps[i]:
                if L[nbp] < nhbest_L[i]:
                    nhbest_L[i] = L[nbp]
                    nhbest_f[i] = f[nbp]
                    nhbest_x[i, :] = x_k[nbp, :]
                    nhbest_i[i] = nbp

    # Initialize stopping criteria distances
    global_dist = 0
    for i in range(swarmsize):
        dist = 0
        for j in range(dimensions):
            dist += (x_k[i, j] - swarm_x[j]) ** 2

        global_dist += dist ** 0.5

    global_distance_reference = global_dist / swarmsize  # relative extent of the swarm

    global_distance = np.zeros(stopIters, float)
    global_L = np.zeros(stopIters, float)
    for k in range(stopIters):
        global_distance[k] = global_distance_reference
        global_L[k] = swarm_L

    # Store History
    if sto_hst:
        logfile.write(rseed, "seed")
        if scale == 1:
            x_uns = np.zeros(x_k.shape)
            for i in range(swarmsize):
                x_uns[i, :] = (x_k[i, :] * space_halflen) + space_centre

        else:
            x_uns = x_k

        if discrete_i:
            for i in range(swarmsize):
                for m in discrete_i:
                    x_uns[i, m] = floor(x_uns[i, m] + 0.5)

        logfile.write(x_uns, "x")
        logfile.write(f, "obj")
        logfile.write(g, "con")
        logfile.write(swarm_x, "gbest_x")
        logfile.write(swarm_f, "gbest_f")
        logfile.write(swarm_g, "gbest_g")

    # Output to Summary File
    if (fileout == 2) or (fileout == 3):
        stext = ""
        stext += "Global Best Particle:\n"
        stext += "-" * 97 + "\n"
        stext += "    Major   Minor   nFCon   Violation(L2)     Objective   Lagrangian   Rel Lagrangian   Global Dist\n"
        stext += "-" * 97 + "\n"
        sfile.write(stext)
        sfile.flush()

    # Outer optimization loop
    k_out = 0
    stop_main_flag = 0
    no_successes = 0
    no_failures = 0
    rho = 1.0
    vcr = 0.0
    while (k_out < maxOutIter) and (stop_main_flag == 0):

        k_out += 1

        # Update g_old Major Iteration
        for i in range(swarmsize):
            g_old[i, :] = g[i, :]

        # Inner optimization loop - core ALPSO algorithm applied to the lagrangian function
        k_inn = 0
        stop_inner = 0
        while (k_inn < maxInnIter) and (stop_inner == 0):

            k_inn += 1

            # calculating new search radius for the best particle ("Guaranteed Convergence" method)
            if (swarm_i == swarm_i_old) and (swarm_L >= swarm_L_old):
                no_failures += 1
                no_successes = 0
            elif (swarm_i == swarm_i_old) and (swarm_L < swarm_L_old):
                no_successes += 1
                no_failures = 0
            else:
                no_successes = 0
                no_failures = 0

            if no_successes > ns:
                rho *= 2.0
                no_successes = 0
            elif no_failures > nf:
                rho *= 0.5
                no_failures = 0

            if rho < 10e-5:
                rho = 10e-5
            elif rho > 1.0:
                rho = 1.0

            # memorization for next outer iteration
            if k_inn == 1:
                swarm_i_old = swarm_i
                swarm_L_old = swarm_L
                swarm_f_old = swarm_f
                swarm_g_old[:] = swarm_g[:]

            # stopping criteria distances
            global_dist = 0
            for i in range(swarmsize):
                dist = 0
                for j in range(dimensions):
                    dist += (x_k[i, j] - swarm_x[j]) ** 2

                global_dist += dist ** 0.5

            global_distance[0] = global_dist / swarmsize  # relative extent of the swarm

            # Update inertia weight
            w = w2 + ((w2 - w1) / global_distance_reference) * global_distance[1]
            if w > w1:
                w = w1
            elif w < w2:
                w = w2

            # Swarm Update
            for i in range(swarmsize):

                # Update velocity vector
                if (nhm == "dlring") or (nhm == "slring") or (nhm == "wheel") or (nhm == "spatial") or (nhm == "sfrac"):
                    lbest_x = nhbest_x[i, :]
                else:
                    lbest_x = swarm_x[:]

                for j in range(dimensions):
                    if i == swarm_i:
                        rr = rand.random()
                        v_k[i, j] = w * v_k[i, j] + -x_k[i, j] + swarm_x[j] + rho * (1.0 - 2.0 * rr)
                    else:
                        r1 = rand.random()
                        r2 = rand.random()
                        rc = rand.random()

                        v_k[i, j] = (
                            w * v_k[i, j]
                            + c1 * r1 * (best_x[i, j] - x_k[i, j]) / dt
                            + c2 * r2 * (lbest_x[j] - x_k[i, j]) / dt
                            + vcr * (1.0 - 2.0 * rc)
                        )

                    # Check for velocity vector out of range
                    if v_k[i, j] > vmax[j]:
                        v_k[i, j] = vmax[j]
                    elif v_k[i, j] < -vmax[j]:
                        v_k[i, j] = -vmax[j]

                    # positions update
                    x_k[i, j] += v_k[i, j] * dt

                    # Check for positions out of range
                    if x_k[i, j] > xmax[j]:
                        x_k[i, j] = xmax[j]
                    elif x_k[i, j] < xmin[j]:
                        x_k[i, j] = xmin[j]

            # Augmented Lagrange
            if h_start:
                [vals, hist_end] = hstfile.read([], ident=["obj", "con"])
                if not hist_end:
                    f = vals["obj"][0]
                    g = vals["con"][0].reshape(g.shape)
                else:
                    h_start = False
                    hstfile.close()

            if not h_start:
                # Evaluate Objective Function
                if scale == 1:
                    xtmp = (x_k * space_halflen) + space_centre
                else:
                    xtmp = x_k

                for m in discrete_i:
                    xtmp[:, m] = floor(xtmp[:, m] + 0.5)

                f, g = objfunc(xtmp)
                nfevals += swarmsize

            # Store History
            if sto_hst:
                if scale == 1:
                    x_uns = np.zeros(x_k.shape)
                    for i in range(swarmsize):
                        x_uns[i, :] = (x_k[i, :] * space_halflen) + space_centre

                else:
                    x_uns = x_k

                if discrete_i:
                    for i in range(swarmsize):
                        for m in discrete_i:
                            x_uns[i, m] = floor(x_uns[i, m] + 0.5)

                logfile.write(x_uns, "x")
                logfile.write(f, "obj")
                logfile.write(g, "con")

            for i in range(swarmsize):

                # Lagrangian Value
                L[i] = f[i]
                if constraints > 0:

                    # Equality Constraints
                    for ell in range(neqcons):
                        tau[i, ell] = g[i, ell]

                    # Inequality Constraints
                    for ell in range(neqcons, constraints):
                        if rp[ell] != 0:
                            if g[i, ell] > -lambda_val[ell] / (2 * rp[ell]):
                                tau[i, ell] = g[i, ell]
                            else:
                                tau[i, ell] = -lambda_val[ell] / (2 * rp[ell])

                        else:
                            tau[i, ell] = g[i, ell]

                    #
                    for ell in range(constraints):
                        L[i] += lambda_val[ell] * tau[i, ell] + rp[ell] * tau[i, ell] ** 2

            # If there is no new better solution for gbest keep the old best position
            # if (L[swarm_i] > swarm_L):
            #    x_k[swarm_i,:] = swarm_x[:]
            #    f[swarm_i] = swarm_f
            #    g[swarm_i,:] = swarm_g[:]
            #    L[swarm_i] = swarm_L

            # Particle Best Update
            for i in range(swarmsize):
                if L[i] < best_L[i]:
                    best_L[i] = L[i]
                    best_f[i] = f[i]
                    best_g[i, :] = g[i, :]
                    best_x[i, :] = x_k[i, :]

            # Swarm Best Update
            for i in range(swarmsize):
                if L[i] < swarm_L:
                    # update of the best particle and best position
                    swarm_i = i
                    swarm_x[:] = x_k[i, :]

                    # update of the best objective function value found
                    swarm_f = f[i]

                    # update of the best constraints values found
                    swarm_g[:] = g[i, :]

                    # update of the swarm best L
                    swarm_L = L[i]

            # Spatial Neighbourhood Update
            if (nhm == "spatial") or (nhm == "sfrac"):
                for i in range(swarmsize):
                    for i2 in range(i + 1, swarmsize):
                        pdist[i, i2] = np.linalg.norm(x_k[i2, :] - x_k[i, :])

                    for i2 in range(i):
                        pdist[i, i2] = pdist[i2, i]

                if nhm == "spatial":
                    for i in range(swarmsize):
                        nhps[i] = []
                        for nb in range(nhn):
                            nhps[i].append(pdist[i, :].argmin())
                            pdist[i, nhps[i][nb]] = inf

                        if nhs == 0:
                            nhps[i].append(i)

                else:
                    frac = ((3 * k_out) + 0.6 * maxOutIter) / maxOutIter
                    if frac >= 1.0:
                        nhm = "gbest"
                    else:
                        for i in range(swarmsize):
                            nhps[i] = []
                            d_max[i] = pdist[i, :].max()
                            for i2 in range(swarmsize):
                                if i == i2:
                                    if nhs == 1:
                                        pass
                                    else:
                                        nhps[i].append(i)

                                else:
                                    if pdist[i, i2] / d_max[i] < frac:
                                        nhps[i].append(i2)

            # Neighbourhood Best Update
            if (nhm == "dlring") or (nhm == "slring") or (nhm == "wheel") or (nhm == "spatial") or (nhm == "sfrac"):
                for i in range(swarmsize):
                    for nbp in nhps[i]:
                        if L[nbp] < nhbest_L[i]:
                            nhbest_L[i] = L[nbp]
                            nhbest_f[i] = f[nbp]
                            nhbest_x[i, :] = x_k[nbp, :]
                            nhbest_i[i] = nbp

            # Print Inner
            if prtInnIter != 0 and np.mod(k_inn, prtInnIter) == 0:
                # output to screen
                print("Outer Iteration: %d     [%d. Inner Iteration]" % (k_out, k_inn))

            if (fileout == 1) or (fileout == 3):
                # output to filename
                pass

            # Inner Loop Convergence
            if k_inn >= minInnIter:
                if swarm_L < swarm_L_old:
                    stop_inner = 1

            # Store History
            if sto_hst:
                logfile.write(swarm_x, "gbest_x")
                logfile.write(swarm_f, "gbest_f")
                logfile.write(swarm_g, "gbest_g")

        # Print Outer
        if prtOutIter != 0 and np.mod(k_out, prtOutIter) == 0:
            # Output to screen
            print("=" * 80 + "\n")
            print("NUMBER OF ITERATIONS: %d\n" % k_out)
            print("NUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" % nfevals)
            print("OBJECTIVE FUNCTION VALUE:")
            print("\tF = %.16g\n" % (float(swarm_f)))
            if constraints > 0:
                # Equality Constraints
                print("EQUALITY CONSTRAINTS VALUES:")
                for ell in range(neqcons):
                    print("\tH(%d) = %g" % (ell, swarm_g[ell]))

                # Inequality Constraints
                print("\nINEQUALITY CONSTRAINTS VALUES:")
                for ell in range(neqcons, constraints):
                    print("\tG(%d) = %g" % (ell, swarm_g[ell]))

            print("\nLAGRANGIAN MULTIPLIERS VALUES:")
            for ell in range(constraints):
                print("\tL(%d) = %g" % (ell, lambda_val[ell]))

            print("\nBEST POSITION:")
            if scale == 1:
                xtmp = (swarm_x[:] * space_halflen) + space_centre
            else:
                xtmp = swarm_x[:]

            for m in discrete_i:
                xtmp[m] = floor(xtmp[m] + 0.5)

            text = ""
            for j in range(dimensions):
                text += "\tP(%d) = %.16g\t" % (j, xtmp[j])
                if np.mod(j + 1, 3) == 0:
                    text += "\n"

            print(text)
            print("=" * 80 + "\n")

        if (fileout == 1) or (fileout == 3):
            # Output to Print File
            ofile.write("\n" + "=" * 80 + "\n")
            ofile.write("\nNUMBER OF ITERATIONS: %d\n" % k_out)
            ofile.write("\nNUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" % nfevals)
            ofile.write("\nOBJECTIVE FUNCTION VALUE:\n")
            ofile.write("\tF = %.16g\n" % (float(swarm_f)))
            if constraints > 0:
                # Equality Constraints
                ofile.write("\nEQUALITY CONSTRAINTS VALUES:\n")
                for ell in range(neqcons):
                    ofile.write("\tH(%d) = %.16g\n" % (ell, swarm_g[ell]))

                # Inequality Constraints
                ofile.write("\nINEQUALITY CONSTRAINTS VALUES:\n")
                for ell in range(neqcons, constraints):
                    ofile.write("\tG(%d) = %.16g\n" % (ell, swarm_g[ell]))

            ofile.write("\nLAGRANGIAN MULTIPLIERS VALUES:\n")
            for ell in range(constraints):
                ofile.write("\tL(%d) = %.16g\n" % (ell, lambda_val[ell]))

            ofile.write("\nPENALTY FACTOR:\n")
            for ell in range(constraints):
                ofile.write("\trp(%d) = %.16g\n" % (ell, rp[ell]))

            ofile.write("\nBEST POSITION:\n")
            if scale == 1:
                xtmp = (swarm_x[:] * space_halflen) + space_centre
            else:
                xtmp = swarm_x[:]

            for m in discrete_i:
                xtmp[m] = floor(xtmp[m] + 0.5)

            text = ""
            for j in range(dimensions):
                text += "\tP(%d) = %.16g\t" % (j, xtmp[j])
                if np.mod(j + 1, 3) == 0:
                    text += "\n"

            ofile.write(text)
            ofile.write("\n" + "=" * 80 + "\n")
            ofile.flush()

        # Store History
        if sto_hst and (minInnIter != maxInnIter):
            logfile.write(k_inn, "ninner")

        # Test Constraint convergence
        stop_con_num = 0
        infeas_con = []
        if constraints == 0:
            stop_constraints_flag = 1
        else:
            for ell in range(neqcons):
                if abs(swarm_g[ell]) <= etol:
                    stop_con_num += 1
                else:
                    infeas_con.append(ell)

            for ell in range(neqcons, constraints):
                if swarm_g[ell] < itol:
                    stop_con_num += 1
                else:
                    infeas_con.append(ell)

            if stop_con_num == constraints:
                stop_constraints_flag = 1
            else:
                stop_constraints_flag = 0

            #        # Test Lagrange multiplier convergence
            #        stop_lambda_flag = 0
            #        if (constraints == 0):
            #            stop_lambda_flag = 1
            #        else:
            #            for l in range(constraints):
            #                if (abs(lambda_val[l]-lambda_old[l]) <= ltol):
            #                    stop_lambda_flag += 1
            #
            #
            #            if (stop_lambda_flag==constraints):
            #                stop_lambda_flag = 1
            #            else:
            #                stop_lambda_flag = 0
            #
            #
            #

        # Test Position and Function convergence
        stop_criteria_flag = 0
        if stopCriteria == 1:

            # setting up the stopping criteria based on distance and tolerance
            for k in range(stopIters - 1, 0, -1):
                global_distance[k] = global_distance[k - 1]
                global_L[k] = global_L[k - 1]

            #
            global_dist = 0
            for i in range(swarmsize):
                dist = 0
                for j in range(dimensions):
                    dist += (x_k[i, j] - swarm_x[j]) ** 2

                global_dist += dist ** 0.5

            global_distance[0] = global_dist / swarmsize  # relative extent of the swarm

            #
            global_L[0] = swarm_L

            #
            if (
                abs(global_distance[0] - global_distance[stopIters - 1]) <= dtol * abs(global_distance[stopIters - 1])
                and abs(global_L[0] - global_L[stopIters - 1]) <= rtol * abs(global_L[stopIters - 1])
                or abs(global_L[0] - global_L[stopIters - 1]) <= atol
            ):
                stop_criteria_flag = 1
            else:
                stop_criteria_flag = 0

        # Test Convergence
        if stop_constraints_flag == 1 and stop_criteria_flag == 1:  # and stop_lambda_flag == 1
            stop_main_flag = 1
        else:
            stop_main_flag = 0

        # Output to Summary File
        if (fileout == 2) or (fileout == 3):
            cvss = 0.0
            for ell in infeas_con:
                cvss += swarm_g[ell] ** 2

            cvL2 = cvss ** 0.5
            if stopCriteria == 1:
                relL = abs(global_L[0] - global_L[stopIters - 1]) / abs(global_L[stopIters - 1])
                stext = "%9d%8d%8d%15.4e%15f%13.4e%16.4e%14.4e\n" % (
                    k_out,
                    k_inn,
                    stop_con_num,
                    cvL2,
                    swarm_f,
                    swarm_L,
                    relL,
                    global_distance[0],
                )
            else:
                stext = "%9d%8d%8d%15.4e%15f%13.4e%16s%14s\n" % (
                    k_out,
                    k_inn,
                    stop_con_num,
                    cvL2,
                    swarm_f,
                    swarm_L,
                    "NA",
                    "NA",
                )

            sfile.write(stext)
            sfile.flush()

        # Update Augmented Lagrangian Terms
        if stop_main_flag == 0:

            if constraints > 0:
                # Update new Tau
                for ell in range(neqcons):
                    tau_new[ell] = swarm_g[ell]

                for ell in range(neqcons, constraints):
                    if swarm_g[ell] > -lambda_val[ell] / (2 * rp[ell]):
                        tau_new[ell] = swarm_g[ell]
                    else:
                        tau_new[ell] = -lambda_val[ell] / (2 * rp[ell])

                # Update Lagrange Multiplier
                for ell in range(constraints):
                    lambda_old[ell] = lambda_val[ell]
                    lambda_val[ell] += 2 * rp[ell] * tau_new[ell]
                    if abs(lambda_val[ell]) < eps:
                        lambda_val[ell] = 0.0

                # Update Penalty Factor
                for ell in range(neqcons):
                    if abs(swarm_g[ell]) > abs(swarm_g_old[ell]) and abs(swarm_g[ell]) > etol:
                        rp[ell] *= 2.0
                    elif abs(swarm_g[ell]) <= etol:
                        rp[ell] *= 0.5

                for ell in range(neqcons, constraints):
                    if swarm_g[ell] > swarm_g_old[ell] and swarm_g[ell] > itol:
                        rp[ell] *= 2.0
                    elif swarm_g[ell] <= itol:
                        rp[ell] *= 0.5

                # Apply Lower Bounds on rp
                for ell in range(neqcons):
                    if rp[ell] < 0.5 * (abs(lambda_val[ell]) / etol) ** 0.5:
                        rp[ell] = 0.5 * (abs(lambda_val[ell]) / etol) ** 0.5

                for ell in range(neqcons, constraints):
                    if rp[ell] < 0.5 * (abs(lambda_val[ell]) / itol) ** 0.5:
                        rp[ell] = 0.5 * (abs(lambda_val[ell]) / itol) ** 0.5

                for ell in range(constraints):
                    if rp[ell] < 1:
                        rp[ell] = 1

            for i in range(swarmsize):
                if constraints > 0:
                    # Update Tau
                    for ell in range(neqcons):
                        tau[i, ell] = g[i, ell]

                    for ell in range(neqcons, constraints):
                        if g[i, ell] > -lambda_val[ell] / (2 * rp[ell]):
                            tau[i, ell] = g[i, ell]
                        else:
                            tau[i, ell] = -lambda_val[ell] / (2 * rp[ell])

            # set craziness velocity for next inner loop run
            vcr = (1 - k_out / maxOutIter) * vcrazy

            # update swarm with new Lagrangian function for next inner run
            for i in range(swarmsize):
                L[i] = f[i]
                if constraints > 0:
                    for ell in range(constraints):
                        L[i] += lambda_val[ell] * tau[i, ell] + rp[ell] * tau[i, ell] ** 2

            swarm_L = L[swarm_i]

            swarm_L_old = swarm_f_old
            if constraints > 0:

                # Equality Constraints
                for ell in range(neqcons):
                    tau_old[ell] = swarm_g_old[ell]

                # Inequality Constraints
                for ell in range(neqcons, constraints):
                    if rp[ell] != 0:
                        if swarm_g_old[ell] > -lambda_val[ell] / (2 * rp[ell]):
                            tau_old[ell] = swarm_g_old[ell]
                        else:
                            tau_old[ell] = -lambda_val[ell] / (2 * rp[ell])

                    else:
                        tau_old[ell] = swarm_g_old[ell]

                #
                for ell in range(constraints):
                    swarm_L_old += lambda_val[ell] * tau_old[ell] + rp[ell] * tau_old[ell] ** 2

            # reset swarm memory for next inner run
            for i in range(swarmsize):
                best_L[i] = L[i]
                best_f[i] = f[i]
                best_g[i, :] = g[i, :]
                best_x[i, :] = x_k[i, :]

    # Print Results
    if prtOutIter != 0:
        # Output to screen
        print("=" * 80 + "\n")
        print("RANDOM SEED VALUE: %.8f\n" % rseed)
        print("NUMBER OF ITERATIONS: %d\n" % k_out)
        print("NUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" % nfevals)
        print("OBJECTIVE FUNCTION VALUE:")
        print("\tF = %.16g\n" % (float(swarm_f)))
        if constraints > 0:
            # Equality Constraints
            print("EQUALITY CONSTRAINTS VALUES:")
            for ell in range(neqcons):
                print("\tH(%d) = %g" % (ell, swarm_g[ell]))

            # Inequality Constraints
            print("\nINEQUALITY CONSTRAINTS VALUES:")
            for ell in range(neqcons, constraints):
                print("\tG(%d) = %g" % (ell, swarm_g[ell]))

        print("\nLAGRANGIAN MULTIPLIERS VALUES:")
        for ell in range(constraints):
            print("\tL(%d) = %g" % (ell, float(lambda_val[ell])))

        print("\nBEST POSITION:")
        if scale == 1:
            xtmp = (swarm_x[:] * space_halflen) + space_centre
        else:
            xtmp = swarm_x[:]

        for m in discrete_i:
            xtmp[m] = floor(xtmp[m] + 0.5)

        text = ""
        for j in range(dimensions):
            text += "\tP(%d) = %.16g\t" % (j, xtmp[j])
            if np.mod(j + 1, 3) == 0:
                text += "\n"

        print(text)
        print("=" * 80 + "\n")

    if (fileout == 1) or (fileout == 3):
        ofile.close()

    if (fileout == 2) or (fileout == 3):
        # Output to Summary
        sfile.write("\n\nSolution:")
        sfile.write("\n" + "=" * 94 + "\n")
        sfile.write("\nNUMBER OF ITERATIONS: %d\n" % k_out)
        sfile.write("\nNUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" % nfevals)
        sfile.write("\nOBJECTIVE FUNCTION VALUE:\n")
        sfile.write("\tF = %.16g\n" % (float(swarm_f)))
        if constraints > 0:
            # Equality Constraints
            sfile.write("\nEQUALITY CONSTRAINTS VALUES:\n")
            for ell in range(neqcons):
                sfile.write("\tH(%d) = %.16g\n" % (ell, swarm_g[ell]))

            # Inequality Constraints
            sfile.write("\nINEQUALITY CONSTRAINTS VALUES:\n")
            for ell in range(neqcons, constraints):
                sfile.write("\tG(%d) = %.16g\n" % (ell, swarm_g[ell]))

        sfile.write("\nLAGRANGIAN MULTIPLIERS VALUES:\n")
        for ell in range(constraints):
            sfile.write("\tL(%d) = %.16g\n" % (ell, float(lambda_val[ell])))

        sfile.write("\nPENALTY FACTOR:\n")
        for ell in range(constraints):
            sfile.write("\trp(%d) = %.16g\n" % (ell, rp[ell]))

        sfile.write("\nBEST POSITION:\n")
        if scale == 1:
            xtmp = (swarm_x[:] * space_halflen) + space_centre
        else:
            xtmp = swarm_x[:]

        for m in discrete_i:
            xtmp[m] = floor(xtmp[m] + 0.5)

        text = ""
        for j in range(dimensions):
            text += "\tP(%d) = %.16g\t" % (j, xtmp[j])
            if np.mod(j + 1, 3) == 0:
                text += "\n"

        sfile.write(text)
        sfile.write("\n" + "=" * 94 + "\n")
        sfile.flush()
        sfile.close()

    # Results
    if scale == 1:
        opt_x = (swarm_x * space_halflen) + space_centre
    else:
        opt_x = swarm_x

    for m in discrete_i:
        opt_x[m] = int(floor(opt_x[m] + 0.5))

    opt_f = swarm_f
    opt_g = swarm_g
    opt_lambda = lambda_val[:]

    return opt_x, opt_f, opt_g, opt_lambda, nfevals, "%.8f" % rseed
