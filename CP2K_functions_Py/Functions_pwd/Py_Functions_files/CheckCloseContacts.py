# calculate close contacts
import numpy as np
clear('all')
BaseFldr = 'G:\Imperial\MattProjects\MD_files\Pit_NVT\Clean\'
system = 'CP_Pit2020_2'
Trajectory = 'CP_Pit2020_2.xyz'
ABC = getABCvectors(BaseFldr,system)
xyz,XYZ,Indx,Atoms,AtomList,nAtoms,startConfig,nConfigs,StepNum = ReadAndParsexyz(BaseFldr,system,Trajectory,ABC)
dist = np.zeros((nAtoms,nAtoms))
for i in np.arange(1,nAtoms+1).reshape(-1):
    for j in np.arange(1,nAtoms+1).reshape(-1):
        if j != i:
            dist[i,j] = np.sqrt(((XYZ(1,i,1) - XYZ(1,j,1)) ** 2) + ((XYZ(1,i,2) - XYZ(1,j,2)) ** 2) + ((XYZ(1,i,3) - XYZ(1,j,3)) ** 2))
        else:
            dist[i,j] = 100

np.amin(dist,[],'all')