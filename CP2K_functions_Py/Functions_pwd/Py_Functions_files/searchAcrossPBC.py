import numpy as np
    
def searchAcrossPBC(Vec = None,Dist = None,XYZ = None,IndxA = None,IndxB = None,ABC = None): 
    PBC_mat = np.array([[1,0,0],[- 1,0,0],[0,1,0],[0,- 1,0],[1,1,0],[- 1,1,0],[1,- 1,0],[- 1,- 1,0],[1,0,1],[- 1,0,1],[0,1,1],[0,- 1,1],[1,1,1],[- 1,1,1],[1,- 1,1],[- 1,- 1,1],[1,0,- 1],[- 1,0,- 1],[0,1,- 1],[0,- 1,- 1],[1,1,- 1],[- 1,1,- 1],[1,- 1,- 1],[- 1,- 1,- 1],[0,0,1],[0,0,- 1]])
    # check for neighbours separated by PBC in x and y only - would be
# good to make this a function, the code is a bit busy here
    for jj in np.arange(1,len(PBC_mat)+1).reshape(-1):
        Vec_PBC = XYZ(IndxB,:) - (XYZ(IndxA,:) + np.multiply(PBC_mat(jj,:),ABC))
        Dist_PBC = np.sqrt((Vec_PBC(:,1) ** 2) + (Vec_PBC(:,2) ** 2) + (Vec_PBC(:,3) ** 2))
        # determine the shortest vector between O atom IndxO(j) and every H after each PBC element is applied
        if np.any(Dist_PBC < Dist):
            PBCIndx = find(Dist_PBC < Dist)
            Vec[PBCIndx,:] = Vec_PBC(PBCIndx,:)
        # determine the shortest distance between O atom IndxO(j) and every H after each PBC element is applied
        Dist = np.amin(Dist_PBC,Dist)
    
    return Vec,Dist
    
    return Vec,Dist