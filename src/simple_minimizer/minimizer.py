# Copyright (C) 2018, 2023 TJ Radcliffe
# Licenced under GPL 3.0

from .vertex import CopyVertex, EmptyVertex, SimpleVertex

import copy
from math import sqrt
import random

# cross product
from numpy import cross

mapConvergenceReasons = {-1: "Exceeded iteration limit", 1: "Closest points indistinguishable", 
                                                2: "Met fractional tolerance", 3:"Minimum scale achieved"}

class SimpleMinimum:
    
    def __init__(self, nAxis, lowerVertex, middleVertex, upperVertex):
        self.nAxis = nAxis
        self.lowerVertex = lowerVertex
        self.middleVertex = middleVertex
        self.upperVertex = upperVertex
        self.fDepth = 0.5*(lowerVertex.fValue+upperVertex.fValue)-middleVertex.fValue
    
    def __lt__(self, other):
        return self.fDepth < other.fDepth
        
    def __str__(self):
        return self.__repr()

    def __repr__(self):
        return " ".join(map(str, (self.nAxis, self.middleVertex, self.fDepth)))

""" Simple minimizer class applies a simple sequential axial minimization.
It does a simple bracketing followed by a parabolic interpolation to 
the minimum along each axis in sequence, repeating with reduced scale
until the scale is small.  This is a crude but extremely effective way 
of decoupling parameters during registration.  If you have an objective
function that goes negative, add a constant to it that is sufficient to
get the minimum around +1.
"""
class SimpleMinimizer:
    
    # Constructor sets dimensions of arrays
    def __init__(self, nDimension):
        
        # The dimension of the space we are working in
        self.nDimension = nDimension

        # The thing we are minimizing
        self.pObjective = None

        # All the points we have sampled
        self.lstHistory = [EmptyVertex(nDimension)]

        # The scale along each axis (used to determine sampling sphere)
        self.lstScales = [1.0 for nI in range(nDimension)]

        # Set the order or lock out dimensions
        self.lstOrder = [nI for nI in range(nDimension)]

        # The points that did not improve things
        self.lstBadPoints = []

        # The fractional radial scale (sampling sphere is this fraction of scale along each axis)
        self.fRadialScale = 1.0

        # The radial scale at which to stop minimizing
        self.fMinimumScale = 0.001
        
        # Fractional tolerance to quit
        self.fFractionalTolerance = 0.001

        # The maximum number of iterations of the outermost loop before giving up
        self.nMaxIterations = 1000

        # The best value of the objective 
        self.fBest = 1e8

        # The second best value
        self.fSecondBest = 1e9

        # The name of the tracefile (cleared after each run, ignored if not set)
        self.strTraceFile = ""

        # The trace file (may be null)
        self.pOutStream = None

    # Reseed the RNG and clear the history vector
    def reseed(self):
      self.lstHistory = []
      self.lstHistory.append(EmptyVertex(self.nDimension))
      random.seed()

    # Objective is what we are minimizing
    def setObjective(self, pObjective):
        self.pObjective = pObjective

    # Set the scale for a given axis.  This determines the initial range of search
    def setScale(self, nAxis, fScale):
        self.lstScales[nAxis] = fScale

    # Disable minimization on a given dimension
    def disableAxis(self, nAxis):
        self.lstMask[nAxis] = -1

    # Enable minimization on a given dimension (will place it in axis order)
    def enableAxis(self, nAxis):
        self.lstMask[nAxis] = nAxis
        
    # Determines when we are close enough to minumum (default is 0.001)
    def setFractionalTolerance(self, fFractionalTolerance):
        self.fFractionalTolerance = fFractionalTolerance
        
    # Set the scales (initial range of search) for all axes
    def setScales(self, lstScales):
        self.lstScales = lstScales

    # Set the order of minimization. -1 means omit. Use carefully!
    def setOrder(self, lstOrder):
        self.lstOrder = lstOrder

    # Reset the order of minimization to the default (in axis order)
    def resetOrder(self, lstOrder):
        self.lstOrder = [nI for nI in range(self.nDimension)]

    # Set the starting position on a given axis
    def setStart(self, nAxis, fStart):
        self.lstHistory[0].setVertex(nAxis, fStart)

    # Set the starting position on all axes
    def setStarts(self, lstStart):
        self.lstHistory[0].lstVertex = lstStart
        
    # Set the expected scale of the minimum as a fraction of the axial scales (quits when reached)
    def setMinimumScale(self, fScale):
        self.fMinimumScale = fScale

    # Set the maximum number of iterations before giving up (default is 1000)
    def setMaxIterations(self, nMaxIterations):
        self.nMaxIterations = nMaxIterations

    # Set the filename for trace output.
    def setTraceFilename(self, strFilename):
        self.strTraceFile = strFilename

    # Get the vertices (last is best)
    def getVertices(self):
        return self.lstHistory

    # Get the best vertex
    def getBestVertex(self):
        return self.lstHistory[-1]
    
    # Get the best value of the objective function
    def getValue(self):
        return self.lstHistory[-1].getValue()

    # Get the best vertex
    def getMinimum(self):
        return self.lstHistory[-1].getVertex()

    # Get reason for convergence
    def getConvergenceReason(self, nReason):
        return mapConvergenceReasons[nReason]

    def adaptiveMinimize(self, lstStart = None):
        """This is still an axial minimizer, but it figures out the order of
        the axes to minimize along rather than requiring the user to do it.
        As always, efficiency is not an option."""
        if lstStart:
            self.setStarts(lstStart)
        
        self.fBest = 1.e8;   # not close to any likely value

        # evaluate starting point
        nCount = 0
        self.lstHistory[-1].setValue(self.pObjective(self.lstHistory[-1].getVertex()))

        # generate collection of minima along each axis
        while True:
            self.fSecondBest = self.fBest;    # record old best value
            lstMinima = []
            lstStart = self.lstHistory[-1]
            for nAxis in range(self.nDimension):            
                lower = EmptyVertex(self.nDimension)
                middle = EmptyVertex(self.nDimension)
                CopyVertex(middle, lstStart)
                upper = EmptyVertex(self.nDimension)
                bFound = self.bracket(nAxis,lower,middle,upper)
                if bFound:
                    lstMinima.append(SimpleMinimum(nAxis, lower, middle, upper))
                # skip axis of not able to bracket
            
            lstMinima.sort()
            lstMinima.reverse()
            for pMinimum in lstMinima:
                nAxis = pMinimum.nAxis
                self.fBest = self.parabolic(nAxis, pMinimum.lowerVertex, pMinimum.middleVertex, pMinimum.upperVertex)        

            self.fRadialScale *= 0.3183098861; # ~1/pi (any ~irrational will do)
            nCount += 1
            if self.nMaxIterations < nCount:   # bail out if we are stuck
                nReason = -1
                break;
                
            # Various reasons for quiting:
            if (self.fBest+self.fSecondBest) == 0:
                nReason = 1
                break
            elif (abs(self.fBest-self.fSecondBest)/(self.fBest+self.fSecondBest) <= self.fFractionalTolerance):
                nReason = 2
                break
            elif (self.fRadialScale <= self.fMinimumScale):
                nReason = 3
                break

        return (nCount, self.lstHistory[-1], nReason);

    # Axial minimization
    def minimize(self, lstStart = None):
        """
        Minimize by minimizations along successive axes.  This is the same
        techinque used in the original pseudo-correlation work, and it has the
        advantage of decoupling the rotation from the translation as much as
        possible. This is really intended for "anatomical" cases where there are
        well-defined axes that are largely independent of each other. The ordering
        of axial minimizations can be controlled using lstOrder, and the 
        system still works pretty well for cases where the axes aren't all
        that independent, but for a lot of cases adaptiveMinimize() may
        work better, as it figures out the optimal ordering of axes as it goes.
        """
        if lstStart:
            self.setStarts(lstStart)
        
        self.fBest = 1.e8;   # not close to any likely value

        # evaluate starting point
        nCount = 0
        self.lstHistory[-1].setValue(self.pObjective(self.lstHistory[-1].getVertex()))

        self.traceOpen()        # open output stream if requested
        self.traceOut()         # dump starting point if we have a trace file
        self.fRadialScale = 1.0
        nReason = 0    # reason for quiting
        while True:
            self.fSecondBest = self.fBest;    # record old best value

            for nAxis in self.lstOrder:
                if nAxis >= 0 and nAxis < self.nDimension:
                    self.fBest = self.minimizeOne(nAxis);    # minimize along the given axis

            self.fRadialScale *= 0.3183098861; # ~1/pi (any ~irrational will do)
            nCount += 1
            if self.nMaxIterations < nCount:   # bail out if we are stuck
                nReason = -1
                break;
                
            # Various reasons for quiting:
            if (self.fBest+self.fSecondBest) == 0:
                nReason = 1
                break
            elif (abs(self.fBest-self.fSecondBest)/(self.fBest+self.fSecondBest) <= self.fFractionalTolerance):
                nReason = 2
                break
            elif (self.fRadialScale <= self.fMinimumScale):
                nReason = 3
                break

        return (nCount, self.lstHistory[-1], nReason);

    # utility for single-axis minimization
    def minimizeOne(self, nAxis):
        """
        Do a brent-based minimization along a given axis.  This is very
        simple-minded, and could be improved by golden section, etc.

            \param nAxis the axis to minimize along
        """
        
        # assume current value really is the best
        fValue = self.lstHistory[-1].getValue();

        # Look for two points that give a triple with
        # a low point in the middle
        lower = EmptyVertex(self.nDimension)
        middle = EmptyVertex(self.nDimension)
        CopyVertex(middle, self.lstHistory[-1])
        upper = EmptyVertex(self.nDimension)
        bFound = self.bracket(nAxis,lower,middle,upper)
        
        if (bFound):
            # at this point we have a triple with the current point the lowest
            # and two higher points on either side.  We do a single parabolic
            # step to try to improve things

            fValue = self.parabolic(nAxis, lower, middle, upper)
        
        return fValue
            
    def parabolic(self, nAxis, lower, middle, upper):
        """Take a single parabolic step"""
        
        fA = lower.getVertex()[nAxis]     # do parabolic fit
        fB = middle.getVertex()[nAxis] 
        fC = upper.getVertex()[nAxis]  
        fFa = lower.getValue() 
        fFb = middle.getValue() 
        fFc = upper.getValue() 

        fBA = fB - fA
        fBC = fB - fC
        fTop = fBA*fBA*(fFb-fFc)-fBC*fBC*(fFb-fFa);
        fBot = fBA*(fFb-fFc)-fBC*(fFb-fFa);

        fDistance = 0.0
        if (0 == fBot):  # NO DIFFERENCE--RETURN EARLY
            self.lstHistory.append(self.lstHistory[-1])
            self.traceOut()  # trace progress if requested
            return self.lstHistory[-1].getValue() 
        else:
            fDistance = -0.5*fTop/fBot     # distance from middle

        middle.incrementVertex(nAxis, fDistance)
        fValue = self.pObjective(middle.getVertex())
        if (fValue > self.lstHistory[-1].getValue()): # no improvement
            fValue = self.lstHistory[-1].getValue()
            self.traceOut() # trace progress if requested
        else:
            middle.setValue(fValue) 
            self.lstHistory.append(middle) 
            self.traceOut()  # trace progress if requested

        return fValue

    # Bracket the minimum along a given axis
    def bracket(self, nAxis, lower, middle, upper):
        """
        Bracket the minimum along an axis using golden section algorithm.  If
        no minimum can be found, return false.  The middle point is passed in
        as a starting point, and passed out as the point in the middle of the
        bracket triple.  The braket routine is based on NRC

            \param nAxis the axis to bracket along
            \param lower the lower (-ve along axis) bracket point
            \param middle the inner point (used to pass in the starting point)
            \param upper the upper (+ve along axis) bracket point

            \return true if minimum found, false otherwise
        """
        bFound = False    # flag to tell if we found anything

        # the distance between any two points
        fDistance = self.fRadialScale*self.lstScales[nAxis]

        if (0 == fDistance): # out of precision
            return False   # NOTE EARLY RETURN (saves needless calc below)
        
        # first point may be higher or lower than middle--start looking along +ve axis
        CopyVertex(upper, middle)
        upper.incrementVertex(nAxis, fDistance)
        upper.setValue(self.pObjective(upper.getVertex()))

        # second point depends on first
        CopyVertex(lower, middle)
        if (upper.getValue() < middle.getValue()):    # keep going
            CopyVertex(lower, middle)
            CopyVertex(middle, upper)
            while True:
                upper.incrementVertex(nAxis, fDistance)
                upper.setValue(self.pObjective(upper.getVertex()))
                if (upper.getValue() < middle.getValue()):   # still going down
                    self.lstHistory.append(upper)
                    self.traceOut()
                    CopyVertex(lower, middle)
                    CopyVertex(middle, upper)
                    fDistance *= 1.618 # expand step by golden ratio and keep trying
                
                else:    # going up again, so we have minimum bracketed
                    bFound = True
                    break
        
        elif (upper.getValue() > middle.getValue()):   # go along -ve axis
        
            while True:
                lower.incrementVertex(nAxis, -fDistance)
                lower.setValue(self.pObjective(lower.getVertex()))
                if (lower.getValue() < middle.getValue()):   # found lower point
                    self.lstHistory.append(lower)
                    self.traceOut()
                    CopyVertex(upper, middle)
                    CopyVertex(middle, lower)
                    fDistance *= 1.618 # expand step by golden ratio and keep trying
                
                else:                
                    bFound = True
                    break
                    
        return bFound  # if points equal, bFound will be false

    # Open the trace file if filename given
    def traceOpen(self):
        if len(self.strTraceFile):
            self.pOutStream = open(self.strTraceFile)

    # Write most recent accepted vertex to the trace file if file open
    def traceOut(self):
        if self.pOutStream:
            self.pOutStream.write(str(self.lstHistory[-1].getValue())+" "+str(self.lstHistory[-1].getVertex())+"\n")
        
    # Close the trace file
    def traceClose(self):
        self.pOutStream = None
