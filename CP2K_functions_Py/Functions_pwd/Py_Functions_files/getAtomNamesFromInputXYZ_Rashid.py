# function [Atoms, AtomList, AtomIndx, PosLines, relevant, nConfigs, nAtoms, lines] = getAtomNamesFromInputXYZ_Rashid(BaseFldr, system)
import numpy as np
    
def getAtomNamesFromInputXYZ_Rashid(BaseFldr = None,system = None): 
    fid = open(np.array([BaseFldr,system,'/',system,'.xyz']))
    lines = textscan(fid,'%s','delimiter','\n','whitespace','')
    # lines = textscan(fid,'#s','delimiter','\n');
    fid.close()
    lines = lines[0]
    nAtoms = str2num(lines[0])
    relevant = find(not cellfun(isempty,strfind(lines,'i =')) )
    # relevant =  find(~cellfun(@isempty,strfind(lines, 'Al')));
    nConfigs = len(relevant)
    PosLines = []
    for i in np.arange(relevant(end()) + 1,relevant(end()) + nAtoms+1).reshape(-1):
        PosLines = np.array([[PosLines],[lines[i]]])
    
    
    # Atoms = strtrim(PosLines(:,1:5));
# AtomList = unique(Atoms,'rows');
# for i = 1:size(AtomList,1)
#     AtomIndx.(strtrim(AtomList(i,:))) = find(sum(Atoms == AtomList(i,:),2) == length(AtomList(i,:)));
# end
    
    return lines,PosLines
    return lines,PosLines