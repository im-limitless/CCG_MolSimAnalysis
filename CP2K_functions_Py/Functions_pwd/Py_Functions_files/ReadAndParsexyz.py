import numpy as np
    
def ReadAndParsexyz(Base = None,Fldr = None,Traj = None,ABC = None): 
    # Base = 'G:\Imperial\MattProjects\MD_files\Pit_NVT\CP_Pit1822\';
# Fldr = 'SS1';
# Traj = 'SS1.xyz';
    
    fid = open(np.array([Base,Fldr,'\',Traj]))
    print('Reading xyz data...')
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    nAtoms = str2num(lines[0])
    relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
    nConfigs = len(relevant)
    Step = []
    Atoms = []
    for i in np.arange(1,nAtoms+1).reshape(-1):
        Atoms = np.array([[Atoms],[lines[(relevant(1) + i)](np.arange(1,5+1))]])
    
    AtomList = unique(Atoms,'rows')
    XYZ = np.zeros((nConfigs,nAtoms,3))
    for i in np.arange(1,AtomList.shape[1-1]+1).reshape(-1):
        setattr(Indx,strtrim(AtomList(i,:)),find(np.sum(Atoms == AtomList(i,:), 2-1) == len(AtomList(i,:))))
        setattr(xyz,strtrim(AtomList(i,:)),np.zeros((nConfigs,len(getattr(Indx,(strtrim(AtomList(i,:))))),3)))
    
    startConfig = 1
    for j in np.arange(startConfig,nConfigs+1).reshape(-1):
        print(np.array(['Processing configuration number ',num2str(j),' of ',num2str(nConfigs)]))
        for k in np.arange(1,nAtoms+1).reshape(-1):
            XYZ[j,k,np.arange[1,3+1]] = str2num(lines[relevant(j) + k](np.arange(5,end()+1)))
        XYZ = wrapXYZ(XYZ,ABC)
        for k in np.arange(1,AtomList.shape[1-1]+1).reshape(-1):
            getattr[xyz,[strtrim[AtomList[k,:]]]][j,:,np.arange[1,3+1]] = XYZ(j,getattr(Indx,(strtrim(AtomList(k,:)))),np.arange(1,3+1))
        if nConfigs > 1:
            initialI = lines[relevant(j)]
            EqI = strfind(initialI,'=')
            ComI = strfind(initialI,',')
            posI = np.array([EqI(1),ComI(1)])
            Step = np.array([[Step],[str2num(initialI(np.arange(posI(1) + 1,posI(2) - 1+1)))]])
    
    return xyz,XYZ,Indx,Atoms,AtomList,nAtoms,startConfig,nConfigs,Step
    return xyz,XYZ,Indx,Atoms,AtomList,nAtoms,startConfig,nConfigs,Step