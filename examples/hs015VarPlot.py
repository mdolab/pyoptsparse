import shelve, numpy,sys
db = {}
db['ipopt']= shelve.open('ipopt_hs015_Hist.hst')
db['slsqp']= shelve.open('slsqp_hs015_Hist.hst')
db['snopt']= shelve.open('snopt_hs015_Hist.hst')
db['fsqp']= shelve.open('fsqp_hs015_Hist.hst')
db['conmin']= shelve.open('conmin_hs015_Hist.hst')
db['nlpql']= shelve.open('nlpql_hs015_Hist.hst')

obj = {}
x1 = {}
x2 = {}

for opt in db.keys():
    n = len(db[opt].keys())
    
    obj[opt] = []
    x1[opt] = []
    x2[opt] = []
    for i in xrange(n):
       
        try:
            obj[opt].append(db[opt]['%d'%i]['fobj'])

            x1[opt].append(db[opt]['%d'%i]['x'][0])
            x2[opt].append(db[opt]['%d'%i]['x'][1])
        except:
            pass
        # end
    # end
# end
from pylab import *

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
styleList=['ko-','ro-','bo-','go-','mo-','co-']
counter=0
for opt in db.keys():
    plot(x1[opt],x2[opt],styleList[counter],label='%s'%(opt))
    counter+=1
# end
plot(xupper,yupper,'k')
legend(loc=3)
xlabel('x1')
ylabel('x2')
title('Simple optimizer comparison')
show()
    
