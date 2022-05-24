import numpy as np
    
def LocateStationaryPoints(f = None): 
    dydx = np.zeros((f.shape[1-1],f.shape[2-1]))
    for i in np.arange(1,f.shape[2-1]+1).reshape(-1):
        #     dydx(:,i) = gradient(smooth(f(:,i))); # smoothed function misses
#     sharp minima...
        dydx[:,i] = gradient(f(:,i))
        count = []
        for j in np.arange(1,len(dydx) - 1+1).reshape(-1):
            if dydx(j) < np.logical_and(0,dydx(j + 1)) > 0:
                count = np.array([count,j])
        MinimaIndx[i] = count
    
    return MinimaIndx
    return MinimaIndx