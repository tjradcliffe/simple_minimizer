import math
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import simple_minimizer

nDimension = 3

class Test:
    """Sphere with origin at [0, 1, 2]"""
    def __init__(self):
        self.lstOrigin = [nI for nI in range(0,nDimension)]
            
    def __call__(self, lstVertex):
        fSum = 0.0
        for nI in range(0, nDimension):
            fSum += (self.lstOrigin[nI]-lstVertex[nI])**2
            
        fSum += 1.0
        
        return math.sqrt(fSum/nDimension)

minimizer = simple_minimizer.SimpleMinimizer(nDimension)
test = Test()
minimizer.setObjective(test)

# test bracketing.  We are starting at [0,0,0] and the minimum
# is at [0, 1, 2]

lower = simple_minimizer.EmptyVertex(nDimension)
middle = simple_minimizer.EmptyVertex(nDimension)
upper = simple_minimizer.EmptyVertex(nDimension)
bFound = minimizer.bracket(0,lower,middle,upper)
print("\nBracket:")
print(((lower.getValue() > middle.getValue()) and (upper.getValue() > middle.getValue())))

print("\nAxial:")
(nCount, vertex, nReason) = minimizer.minimize()
lstVertex = vertex.getVertex()
fError = math.sqrt(sum([(nI-lstVertex[nI])**2 for nI in range(0,nDimension)])/nDimension)
print("Error, reason: ", fError, nReason)
print(fError < 0.001)
print(nCount < 10)

print("\nAdaptive:")
minimizer.reseed()
(nCount, vertex, nReason) = minimizer.adaptiveMinimize()
lstVertex = vertex.getVertex()
fError = math.sqrt(sum([(nI-lstVertex[nI])**2 for nI in range(0,nDimension)])/nDimension)
print("Error, reason: ", fError, nReason)
print(fError < 0.001)
print(nCount < 10)
