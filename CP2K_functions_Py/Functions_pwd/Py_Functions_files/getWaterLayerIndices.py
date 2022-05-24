import numpy as np
    
def getWaterLayerIndices(Indx = None,XYZ = None,Dens_O = None,z = None): 
    Minima = LocateStationaryPoints(mean(Dens_O,2))
    for i in np.arange(1,XYZ.shape[1-1]+1).reshape(-1):
        FirstLayerIndx[i] = np.array([[intersect(Indx.O,find(XYZ(i,:,3) <= z(Minima[0](1))))],[intersect(Indx.O,find(XYZ(i,:,3) >= z(Minima[0](end()))))]])
        SecondLayerIndx[i] = np.array([[intersect(Indx.O,find(XYZ(i,:,3) > np.logical_and(z(Minima[0](1)),XYZ(i,:,3)) < z(Minima[0](2))))],[intersect(Indx.O,find(XYZ(i,:,3) < np.logical_and(z(Minima[0](end())),XYZ(i,:,3)) > z(Minima[0](end() - 1))))]])
    
    return FirstLayerIndx,SecondLayerIndx
    return FirstLayerIndx,SecondLayerIndx