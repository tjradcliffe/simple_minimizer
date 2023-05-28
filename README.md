# Simple Minimizer

This minimizer was designed to be brutally simple and extremely
robust. The original aim was primarily multi-model image registration
during the development of the pseudo-correlation image registration
algorithm for online portal imaging and other applications.

It uses a combination of bracketing and parabolic interpolation to
get the job done. It is not designed for speed.

An understanding of the shape of your objective function is
desirable. It is useful to map it on various sets of axes. 

The axes you use to represent your problem are not in general 
independent: your objective will have diagonal "troughs"
due to the ability of one axis to trade off against another.

In general you are not minimizing in a vector space, and if you are the
scales on different axes are often wildly different.

The value of the objective function should be strictly positive over the
domain, and ideally around 1.0 near the minimum, where "around"
means 1E-4 => 1E4 or so.

If local minima are a problem call `.reseed()` on your minimizer object
after recording the current minimum, set new starting points (or use
the one you've just found) and re-run. For very hard problems with a
lot of local minimia I've sometimes run this three or five times and looked
for a majority vote on the true minimum or a lowest value (which works
best depends on the scale of the noise on the objective function: if the
noise is smooth on the scale of the minimum lowest value works well, if
the noise is rough on the scale of the minimum majority of minima within
a small region works well.)

Simple usage example:

```
from math import sqrt

from simple_minimizer import *

nDimension = 3

# objective function is a callable
class Test(object):
    def __init__(self):
        # place we are looking for
        self.lstOrigin = [nI for nI in range(0,nDimension)]
            
    def __call__(self, lstVertex):
        # lstVertex is a list of floats giving location on axes
        fSum = 0.0
        for nI in range(0, nDimension):
            fSum += (self.lstOrigin[nI]-lstVertex[nI])**2
            
        fSum += 1.0 # keeping it positive and ~ 1
        
        return sqrt(fSum/nDimension)

minimizer = SimpleMinimizer(nDimension)
test = Test()
minimizer.setObjective(test)

# nReason = -1 => failed to converge, 1 => best/second-best equal, 2 => ratio < 0.001, 3 => minimum scale
(nCount, vertex, nReason) = minimizer.minimize()
lstVertex = vertex.getVertex() # extract the list of floats from the vertex object
fValue = vertex.getValue() # the value of the objective function at the minimum
fError = sqrt(sum([(nI-lstVertex[nI])**2 for nI in range(0,nDimension)])/nDimension)
print(fError < 0.001)
print(nCount < 10)
print(lstVertex)
    
```

In a realistic cases the objective function will be ferociously complex. Some
of mine have involved computing DRRs (digitially reconstructed radiograms
from CT or MR data) on each call.
