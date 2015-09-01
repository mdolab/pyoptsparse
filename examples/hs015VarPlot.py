import sys
import numpy as np
import matplotlib.pyplot as plt
from pyoptsparse.pyoptsparse.sqlitedict.sqlitedict import SqliteDict

db = SqliteDict()
opts = ['ipopt', 'slsqp', 'snopt', 'fsqp', 'conmin', 'nlpqlp', 'psqp']
for opt in opts:
    fileName = '%s_hs015_Hist.hst'%opt
    try:
        db[opt] = SqliteDict(fileName)
    except:
        pass

obj = {}
x1 = {}
x2 = {}

for opt in db.keys():
    n = int(db[opt]['last'])
    
    obj[opt] = []
    x1[opt] = []
    x2[opt] = []
    for i in xrange(n):
        try:
            obj[opt].append(db[opt]['%d'%i]['funcs']['obj'])
            x1[opt].append(db[opt]['%d'%i]['xuser']['xvars'][0])
            x2[opt].append(db[opt]['%d'%i]['xuser']['xvars'][1])
        except:
            pass
   
# Generate the Rosenbrock contours
delta = 0.25
x = np.arange(-2.5, 1.5, delta)
y = np.arange(-6.5, 3.0, delta)
X, Y = np.meshgrid(x, y)

Z = 100*(Y - X**2)**2 + (1-X)**2
# and the constraint contours
A = X*Y
B = X + Y**2

# plot the contours and constraints
plt.figure(figsize=(11, 8.5))
levels = [6000,3000,2000,1000,500,250]
CS = plt.contour(X, Y, Z,levels,colors='k')
levels = np.arange(1.0, 1.01)
CS1 = plt.contour(X,Y,A, levels,colors='g')
levels = np.arange(0.0, 0.01)
CS2 = plt.contour(X,Y,B, levels,colors='b')
plt.clabel(CS, inline=1, fontsize=10)
plt.clabel(CS1, inline=1, fontsize=10)
plt.clabel(CS2, inline=1, fontsize=10)

# set the one sided variable
xupper=[0.5,0.5]
yupper = [-7,3.0]

# Now plot the optimizer output
styleList=['ko-','ro-','bo-','go-','mo-','co-','ks--']
counter=0
for opt in db.keys():
    plt.plot(x1[opt],x2[opt],styleList[counter],label='%s'%(opt))
    counter+=1

# end
plt.plot(xupper,yupper,'k')
plt.legend(loc=3)
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('Simple optimizer comparison')
plt.show()
