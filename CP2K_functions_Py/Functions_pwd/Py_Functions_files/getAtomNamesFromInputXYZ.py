import numpy as np
    
def getAtomNamesFromInputXYZ(BaseFldr = None,system = None): 
    fid = open(np.array([BaseFldr,system,'\',system,'.xyz']))
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    fid.close()
    lines = lines[0]
    nAtoms = str2num(lines[0])
    relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
    nConfigs = len(relevant)
    PosLines = []
    for i in np.arange(relevant(end()) + 1,relevant(end()) + nAtoms+1).reshape(-1):
        PosLines = np.array([[PosLines],[lines[i]]])
    
    Atoms = strtrim(PosLines(:,np.arange(1,5+1)))
    AtomList = unique(Atoms,'rows')
    for i in np.arange(1,AtomList.shape[1-1]+1).reshape(-1):
        setattr(AtomIndx,strtrim(AtomList(i,:)),find(np.sum(Atoms == AtomList(i,:), 2-1) == len(AtomList(i,:))))
    
    return Atoms,AtomList,AtomIndx
    return Atoms,AtomList,AtomIndx