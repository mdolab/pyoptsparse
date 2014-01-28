# Define all the possible sets of parameters we have.
import os
sens = ['FD','CS','none','user']
sensMode = ['', 'pgc']
constrained =[0, 1]
useDict = [0, 1]
testHist = ['no', 'hot', 'cold']
groups = [0, 1]

nproc = 1
for sm in sensMode:
    for s in sens:
        for c in constrained:
            for d in useDict:
                for h in testHist:
                    for g in groups:
                        if sm == 'pgc':
                            nproc = 2
                        else:
                            nproc = 1 
                        cmd = 'mpirun -np %d python rosenbrock.py --sensMode=%s --sens=%s --constrained=%s --useDict=%s --testHist=%s --groups=%s'%(nproc, sm, s,c,d,h,g)
                            
                        print '-'*120
                        print cmd
                        print '-'*120
                        os.system(cmd)
