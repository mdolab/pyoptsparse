'''
This script parses the output of SNOPT Printout files. It extracts the optimality, feasibility and merit function for each major iteration, and writes a tecplot file that can be used for the generation of plots.

Usage: run it as
python SNOPT_parse.py path_to_file

If no filename is given, the default 'SNOPT_print.out' will be assumed.
'''

import sys
import numpy as np

if len(sys.argv) > 1:
    file_name = sys.argv[1]
else:
    file_name = 'SNOPT_print.out'

# Values we want to extract
merit = []
optimality = []
feasibility = []

# Open file
f = open(file_name,'r')

QP_line = "    Itn       QPmult  QPstep   nInf   SumInf   rgNorm    QPobjective   +SBS   -SBS    -BS    Pivot     L+U ncp    nS  condHz"
LP_line = "    Itn       LPmult  LPstep   nInf   SumInf             LPobjective   +SBS   -SBS    -BS    Pivot     L+U ncp    nS"
NP_line = "   Itns Major Minors    Step   nCon Feasible  Optimal  MeritFunction     L+U BSwap     nS  condHz Penalty"
EXIT_line = "SNOPTC EXIT"

# Read all lines
read_next_line = False
for line in f:

    if EXIT_line in line:
        break

    if QP_line in line or LP_line in line:
        read_next_line = False

    if NP_line in line:
        read_next_line = True
        continue

    if read_next_line and len(line) > 100:
        merit.append(float(line[55:68]))
        optimality.append(float(line[46:53]))
        feasibility.append(float(line[37:44]))

# Close print file
f.close()

# Now dump out the data into a tecplot file:
f = open("SNOPT.dat",'w')
iterations = np.arange(len(merit))
f.write ('VARIABLES = "Iteration","Merit Function","Optimality","Feasibility",\n')
f.write('Zone T="SNOPT Data" I=%d\n'%(len(merit)))
f.write('DATAPACKING=POINT\n')
for i in range(len(merit)):
    fea = max(feasibility[i],1e-6)
    f.write('%f %f %f %f\n'%(iterations[i],merit[i],optimality[i],fea))
f.close()