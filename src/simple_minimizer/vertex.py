
# Copyright (C) 2018 TJ Radcliffe 
# Licenced under GPL 3.0

def EmptyVertex(nDimensions):
    return SimpleVertex([0.0 for nI in range(0, nDimensions)], 0.0, 0.0)

def CopyVertex(v1, v2):
    v1.lstVertex = [v2.lstVertex[nI] for nI in range(0, len(v2.lstVertex))]
    v1.fValue = v2.fValue
    v1.fDistance = v2.fDistance

# Helper class to hold vertex data
class SimpleVertex:
    
    def __init__(self, lstVertex, fValue, fDistance):
        
        # The position of the vertex
        self.lstVertex = [lstVertex[nI] for nI in range(0,len(lstVertex))]

        # The value at the vertex
        self.fValue = fValue

        # The distance between this vertex and the previous one
        self.fDistance = fDistance
        
    def __repr__(self):
        return str(self.fValue)+" "+str(self.lstVertex)
        
    def __getitem__(self, nIndex):
        return self.lstVertex[nIndex]
        
    def getVertex(self):
        return self.lstVertex
        
    def setVertex(self, nAxis, fDistance):
        self.lstVertex[nAxis] = fDistance
        
    def incrementVertex(self, nAxis, fDistance):
        self.lstVertex[nAxis] += fDistance

    def getValue(self):
        return self.fValue

    def setValue(self, fValue):
        self.fValue = fValue

    def getDistance(self):
        return self.fDistance

if __name__ == "__main__":

    nSize = 3
    lstVertex = [nI for nI in range(0, nSize)]
    fDistance = 10.5
    fValue = 1.234
    vertex = SimpleVertex(lstVertex, fValue, fDistance)
    print(nSize == len(vertex.getVertex()))
    print(fValue == vertex.getValue())
    print(fDistance == vertex.getDistance())
    
    nSize = 5
    empty_vertex = EmptyVertex(nSize)
    print(nSize == len(empty_vertex.getVertex()))
    print(0.0 == empty_vertex.getValue())
    print(0.0 == empty_vertex.getDistance())

