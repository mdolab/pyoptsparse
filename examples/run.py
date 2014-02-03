# Define all the possible sets of parameters we have.
import os
sens = ['user','FD','CS']#,'none']
sensMode = ['']#, 'pgc']
constrained =[1, 0]
testHist = ['no']#, 'hot', 'cold']
groups = [0, 1]

nproc = 1
for sm in sensMode:
    for s in sens:
        for c in constrained:
            for h in testHist:
                for g in groups:
                    if sm == 'pgc':
                        nproc = 2
                    else:
                        nproc = 1 
                    cmd = 'mpirun -np %d python rosenbrock.py --sensMode=%s --sens=%s --constrained=%s --testHist=%s --groups=%s'%(nproc, sm, s,c,h,g)
                            
                    print '-'*120
                    print cmd
                    print '-'*120
                    os.system(cmd)
